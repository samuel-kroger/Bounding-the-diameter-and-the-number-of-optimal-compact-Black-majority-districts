"""
run_gerrychain.py -- heuristic UPPER BOUNDS on the diameter (s) vs. #majority-
minority tradeoff, via GerryChain short-bursts + diameter minimization.

For each instance we:
  Phase 1 (find MM): Gingleator short-bursts (the cited SOTA for maximizing the
    number of majority-minority districts) -- discover the largest #MM reachable
    and keep, for each #MM level, the smallest-diameter plan seen.
  Phase 2 (shrink diameter): for each #MM target t = 1..max, run a
    SingleMetricOptimizer that MINIMIZES the max within-district graph diameter,
    seeded from a phase-1 MM plan and with an MM-floor penalty so it cannot drift
    below t. The best diameter is an UPPER BOUND on s*(t).

Any plan returned has >= t MM districts and a known diameter d, so:
  * s*(t) <= d   (a heuristic upper bound, complementing the paper's Section-6
                  lower bounds -- together they bracket the tradeoff), and
  * the plan is a feasible WARM START for the exact model at s = d.

Outputs (into ../Polsby_Popper_optimization-main/results_gerrychain_<date>/ so
score_belotti_plans.py picks them up):
  * <state>_<level>_gerrychain_ub.csv         -- best plan (max #MM, min diam) as
                                                 GEOID20,district, for scoring +
                                                 warm starting.
  * gerrychain_tradeoff.csv                    -- (state, level, #MM, min diameter)
                                                 rows: the upper-bound curve.

Population deviation is matched to the exact model (+/-0.5%, i.e. 1% total), so
the plans are valid in the same feasible region.

Requires: pip install gerrychain  (latest; 0.3.x).
Run from code/:  python run_gerrychain.py
"""
import csv
import json
import os
import sys
import time
from datetime import date
from functools import partial

import networkx as nx
import fast_graph as fg

try:
    from gerrychain import Partition, Graph, constraints, proposals
    from gerrychain.optimization import Gingleator, SingleMetricOptimizer
    from gerrychain.tree import recursive_tree_part, bipartition_tree
    from gerrychain.updaters import Tally, cut_edges
except Exception as e:  # pragma: no cover
    print("GerryChain is required: pip install gerrychain  (latest).  (%s)" % e)
    sys.exit(1)

# Seed bipartition: just many attempts (NO pair reselection -- recursive_tree_part
# cannot handle the ReselectException it raises). High max_attempts lets tight
# (+/-0.5%) deviation seed on small county graphs where the default 10k fails.
SEED_BIPARTITION = partial(bipartition_tree, max_attempts=100000)
# ReCom proposal bipartition: pair reselection IS appropriate here (the proposal
# catches ReselectException and picks a different district pair).
ROBUST_BIPARTITION = partial(bipartition_tree, allow_pair_reselection=True,
                             max_attempts=100000)

# (state, level, plural, group, k)
INSTANCES = [
    ("MS", "county", "counties", "black", 4),
    ("AL", "county", "counties", "black", 7),
    ("LA", "county", "counties", "black", 6),
    ("SC", "county", "counties", "black", 7),
    ("NM", "tract", "tracts", "hispanic", 3),
    ("MS", "tract", "tracts", "black", 4),
    ("SC", "tract", "tracts", "black", 7),
    ("LA", "tract", "tracts", "black", 6),
    ("CO", "tract", "tracts", "hispanic", 8),
]

DEVIATION = 0.01          # total; +/- DEVIATION/2, matches the exact model
EPS = DEVIATION / 2.0     # GerryChain "within_percent_of_ideal_population"
NODE_REPEATS = 20         # ReCom cut-finding persistence (higher = more robust)
PHASE1_BURST_LEN = 10
PHASE1_NUM_BURSTS = 300   # phase-1 plans = burst_len * num_bursts (find MM)
PHASE2_BURST_LEN = 10
PHASE2_NUM_BURSTS = 200   # phase-2 plans per #MM target (shrink diameter)
# Override both burst counts from the command line with  bursts=N  (see main()).

BLACK_ANYPART = [
    "P0030004", "P0030011", "P0030016", "P0030017", "P0030018", "P0030019",
    "P0030027", "P0030028", "P0030029", "P0030030", "P0030037", "P0030038",
    "P0030039", "P0030040", "P0030041", "P0030042", "P0030048", "P0030049",
    "P0030050", "P0030051", "P0030052", "P0030053", "P0030058", "P0030059",
    "P0030060", "P0030061", "P0030064", "P0030065", "P0030066", "P0030067",
    "P0030069", "P0030071",
]


def load_graph(state, level, plural, group):
    with open("raw_data/%s/json/%s_%s.json" % (level, state, plural)) as f:
        d = json.load(f)
    G = nx.Graph()
    geoid = {}
    for nd in d["nodes"]:
        i = nd["id"]
        mvap = sum(nd[c] for c in BLACK_ANYPART) if group == "black" else nd["P0040002"]
        G.add_node(i, TOTPOP=nd["P0010001"], VAP=nd["P0030001"], MVAP=mvap)
        geoid[i] = str(nd["GEOID20"])
    for i, nbrs in enumerate(d["adjacency"]):
        for e in nbrs:
            G.add_edge(i, e["id"])
    return Graph(G), geoid


def make_diam_fn(G, csr):
    n = G.number_of_nodes()
    member = [False] * n

    def maxdiam(parts):
        worst = 0
        for key in parts:
            D = list(parts[key])
            for v in D:
                member[v] = True
            for a in D:
                g = csr.bfs1(a, allowed=member)
                for b in D:
                    if csr._s1[b] == g and csr._d1[b] > worst:
                        worst = csr._d1[b]
            for v in D:
                member[v] = False
        return worst

    return maxdiam


def n_mm(part):
    return sum(1 for d in part.parts if part["MVAP"][d] >= 0.5 * part["VAP"][d])


def run_instance(state, level, plural, group, k, outdir, tradeoff_rows):
    print("\n=== %s %s (%s), k=%d ===" % (state, level, group, k))
    G, geoid = load_graph(state, level, plural, group)
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    maxdiam = make_diam_fn(G, csr)
    ideal = sum(G.nodes[i]["TOTPOP"] for i in G.nodes) / k

    upd = {"population": Tally("TOTPOP", alias="population"),
           "MVAP": Tally("MVAP"), "VAP": Tally("VAP"), "cut_edges": cut_edges}

    # Build the seed partition with a persistent bipartition; if +/-EPS is too
    # tight for ReCom to even seed (common for small county graphs), escalate the
    # deviation and note it. Matching the exact model's +/-0.5% is preferred, but
    # a feasible heuristic plan at a looser deviation is better than none.
    assignment, use_eps = None, EPS
    for mult in (1, 2, 4, 8):
        e = EPS * mult
        try:
            assignment = recursive_tree_part(G, range(k), ideal, "TOTPOP", e,
                                             NODE_REPEATS, method=SEED_BIPARTITION)
            use_eps = e
            break
        except Exception:
            continue
    if assignment is None:
        print("  skipped: ReCom could not seed a balanced plan even at +/-%.1f%% "
              "deviation (instance too tight for the chain)." % (100 * EPS * 8))
        return
    if use_eps > EPS:
        print("  NOTE: using +/-%.2f%% deviation here (the exact model's +/-%.2f%% "
              "is too tight for ReCom to seed on this instance)."
              % (100 * use_eps, 100 * EPS))
    init = Partition(G, assignment, upd)
    prop = partial(proposals.recom, pop_col="TOTPOP", pop_target=ideal,
                   epsilon=use_eps, node_repeats=NODE_REPEATS, method=ROBUST_BIPARTITION)
    cc = [constraints.within_percent_of_ideal_population(init, use_eps), constraints.contiguous]

    t0 = time.time()

    # ---- Phase 1: maximize #MM, keep min-diameter plan per #MM level ----
    best_diam = {}      # mm -> min diameter
    best_assign = {}    # mm -> assignment dict
    gg = Gingleator(prop, cc, init, minority_pop_col="MVAP", total_pop_col="VAP",
                    threshold=0.5, score_function=Gingleator.reward_partial_dist)
    for part in gg.short_bursts(PHASE1_BURST_LEN, PHASE1_NUM_BURSTS, with_progress_bar=False):
        mm = n_mm(part)
        d = maxdiam(part.parts)
        if mm not in best_diam or d < best_diam[mm]:
            best_diam[mm] = d
            best_assign[mm] = dict(part.assignment)
    max_mm = max(best_diam) if best_diam else 0
    print("  phase 1: reached up to #MM=%d (diam %s)"
          % (max_mm, best_diam.get(max_mm)))

    # ---- Phase 2: minimize diameter for each #MM target, seeded ----
    for target in range(1, max_mm + 1):
        # seed from the best phase-1 plan that already has >= target MM
        seed_mm = min((m for m in best_assign if m >= target), default=None)
        if seed_mm is None:
            continue
        seed = Partition(G, best_assign[seed_mm], upd)

        def obj(p, _t=target):
            return maxdiam(p.parts) + 100000 * max(0, _t - n_mm(p))

        opt = SingleMetricOptimizer(prop, cc, seed, optimization_metric=obj, maximize=False)
        for part in opt.short_bursts(PHASE2_BURST_LEN, PHASE2_NUM_BURSTS, with_progress_bar=False):
            if n_mm(part) >= target:
                d = maxdiam(part.parts)
                if d < best_diam.get(target, 1 << 30):
                    best_diam[target] = d
                    best_assign[target] = dict(part.assignment)
        print("  phase 2: #MM>=%d -> min diameter (s* upper bound) = %d"
              % (target, best_diam[target]))
        tradeoff_rows.append(dict(state=state, level=level, mm=target,
                                  min_diameter=best_diam[target]))

    # ---- save the representative plan (max #MM, its min diameter) ----
    if max_mm >= 1:
        rep = max_mm
        assign = best_assign[rep]
        path = os.path.join(outdir, "%s_%s_gerrychain_ub.csv" % (state, level))
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["GEOID20", "district"])
            for i in range(n):
                w.writerow([geoid[i], assign[i]])
        print("  saved best plan (#MM=%d, diam=%d) -> %s  (%.1fs)"
              % (rep, best_diam[rep], path, time.time() - t0))
    else:
        print("  no MM plan found (%.1fs)" % (time.time() - t0))


def main():
    """CLI:
        python run_gerrychain.py                      # all instances, default bursts
        python run_gerrychain.py SC:tract LA:tract    # only these instances
        python run_gerrychain.py SC:tract bursts=600  # these, with 600 bursts/phase
        python run_gerrychain.py bursts=600           # all, 600 bursts/phase
    Re-running a subset MERGES into the existing tradeoff file (other instances'
    rows and saved plans are preserved), so you can push the stubborn tracts
    harder without redoing the easy ones."""
    global PHASE1_NUM_BURSTS, PHASE2_NUM_BURSTS
    wanted = []
    for a in sys.argv[1:]:
        if a.startswith("bursts="):
            PHASE1_NUM_BURSTS = PHASE2_NUM_BURSTS = int(a.split("=", 1)[1])
        else:
            wanted.append(a.replace("_", ":").upper())

    instances = INSTANCES
    if wanted:
        instances = [t for t in INSTANCES if ("%s:%s" % (t[0], t[1])).upper() in wanted]
        if not instances:
            print("No instances matched %s.  Use e.g.:  SC:tract LA:tract [bursts=600]"
                  % " ".join(sys.argv[1:]))
            return
    print("Running %d instance(s), %d/%d bursts (phase1/phase2): %s"
          % (len(instances), PHASE1_NUM_BURSTS, PHASE2_NUM_BURSTS,
             ", ".join("%s %s" % (t[0], t[1]) for t in instances)))

    today = date.today().strftime("%Y_%b_%d")
    outdir = os.path.join("..", "Polsby_Popper_optimization-main",
                          "results_gerrychain_" + today)
    os.makedirs(outdir, exist_ok=True)

    new_rows = []
    for (state, level, plural, group, k) in instances:
        try:
            run_instance(state, level, plural, group, k, outdir, new_rows)
        except Exception as e:
            import traceback
            print("  ERROR on %s %s: %s" % (state, level, e))
            traceback.print_exc()

    # merge into gerrychain_tradeoff.csv: keep rows for instances we did NOT run,
    # overwrite (state, level, mm) rows for the ones we did.
    tpath = os.path.join(outdir, "gerrychain_tradeoff.csv")
    ran = {(t[0], t[1]) for t in instances}
    merged = {}
    if os.path.exists(tpath):
        for r in csv.DictReader(open(tpath)):
            if (r["state"], r["level"]) not in ran:
                merged[(r["state"], r["level"], str(r["mm"]))] = r
    for r in new_rows:
        merged[(r["state"], r["level"], str(r["mm"]))] = r
    with open(tpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["state", "level", "mm", "min_diameter"])
        w.writeheader()
        for r in merged.values():
            w.writerow(r)
    print("\nWrote plans + gerrychain_tradeoff.csv to", outdir)
    print("Now run:  python run_A1_ours.py   then   python score_belotti_plans.py")


if __name__ == "__main__":
    main()
