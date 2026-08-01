"""
run_A1_ours.py  --  run OUR side of the A1 benchmark (full optimized diameter
model: reduced formulation + parcel fixing + independent-set symmetry-breaking +
clique distance cuts + z-ordering).

For each instance we now climb s until the model yields at least TARGET_MM
majority-minority districts (default 1), rather than stopping at the first
*feasible* plan (which can have #MM = 0).  Concretely, at each s we solve the
MM-maximizing model with diameter <= s and:

  * #MM (optimal) >= target              -> success, stop.
  * #MM (optimal)  < target              -> larger s can only help, so s += 1.
  * Infeasible (no valid districting)    -> s += 1.
  * Timeout/other with an incumbent      -> if the incumbent already has
                                            #MM >= target, stop; else s += 1.
  * Timeout/other with NO incumbent      -> s += 1 (larger s is easier to be
                                            feasible), logging the reason.

The climb is capped by the graph's upper bound on s and by a cumulative
INSTANCE_TIME_BUDGET so a hopeless instance cannot loop forever.  EVERY instance
emits exactly one summary row -- including ones that captured nothing -- with a
human-readable reason in the "compact" column, and full tracebacks are printed
for any exception (so the SC-county / tract / LA-tract failures explain
themselves on the re-run).

Run from the code/ directory (needs gurobipy + the project pipeline):
    python run_A1_ours.py
"""
import csv
import glob
import os
import time
import traceback
from collections import defaultdict
import fast_graph as fg
import benchmark_run_helpers as H
from benchmark_to_excel import build_workbook
from comp_experiment import problem_instance


def load_gerrychain_warmstart(state, level, G, mvap, vap, f=0.5):
    """Find a heuristic plan for this instance and turn it into a warm start:
    ((minority_districts, majority_districts), diam), with majority-minority
    districts (any-part MVAP >= f*VAP) listed FIRST and the plan's max
    within-district graph diameter. Prefers the GerryChain plan (tightest
    diameter); falls back to the Belotti CCC plan when GerryChain found no MM plan
    (e.g. LA tract). Returns None if no usable plan exists."""
    cands = glob.glob(os.path.join(
        "..", "Polsby_Popper_optimization-main",
        "results_gerrychain_*", "%s_%s_gerrychain_ub.csv" % (state, level)))
    src = "GerryChain"
    if not [c for c in cands if os.path.exists(c)]:
        # fall back to the CCC plan (also a valid +/-0.5% MM plan)
        cands = glob.glob(os.path.join(
            "..", "Polsby_Popper_optimization-main",
            "results_*", "%s_%s_gingles_ccc.csv" % (state, level)))
        src = "CCC"
    cands = [c for c in cands if os.path.exists(c)]
    if not cands:
        return None
    path = max(cands, key=os.path.getmtime)
    g2i = {}
    for v in G.nodes():
        gid = G.nodes[v].get("GEOID20")
        if gid is not None:
            g2i[str(gid).strip()] = v
    n = G.number_of_nodes()
    labels = [None] * n
    for r in csv.DictReader(open(path)):
        idx = g2i.get(str(r["GEOID20"]).strip())
        if idx is not None:
            labels[idx] = int(r["district"])
    if any(l is None for l in labels):
        print("  (GerryChain plan %s did not match all nodes; skipping warm start)" % path)
        return None
    dd = defaultdict(list)
    for v, j in enumerate(labels):
        dd[j].append(v)
    districts = list(dd.values())
    mm_dists, non_dists = [], []
    for D in districts:
        m = sum(mvap[v] for v in D)
        vp = sum(vap[v] for v in D)
        (mm_dists if (vp and m >= f * vp) else non_dists).append(D)
    diam = H.plan_diameter(G, districts)
    print("  warm start (%s): %s  #MM=%d  diameter=%d" % (src, path, len(mm_dists), diam))
    return (mm_dists, non_dists), diam

# (state, level, group, target_mm)  -- target_mm is the smallest #MM we want a
# returned plan to have (so the row is a real majority-minority comparison).
INSTANCES = [
    ("MS", "county", "black",    1),
    ("AL", "county", "black",    1),
    ("LA", "county", "black",    1),
    ("SC", "county", "black",    1),
    ("NM", "tract",  "hispanic", 1),
    ("MS", "tract",  "black",    1),
    ("SC", "tract",  "black",    1),
    ("LA", "tract",  "black",    1),
    ("CO", "tract",  "hispanic", 1),
]

PER_MODEL_TIME_LIMIT = 3600        # one hour per MIP solve (per s)
INSTANCE_TIME_BUDGET = 8 * 3600    # stop climbing s once cumulative solve time
                                   # exceeds this (prevents runaway on hopeless
                                   # instances). Raised to 8h so the tract climbs
                                   # reach the (large) s where an MM district fits.


def _gap(m):
    try:
        return m.MIPGap
    except Exception:
        return None


def _objbound(m):
    try:
        return m.ObjBound
    except Exception:
        return None


def _err_row(state, level, n, k, note):
    return dict(state=state, level=level, n=n, k=k,
                method="diameter (ours, full)", compact=note,
                mm=None, time=None, gap=None,
                pp=None, reock=None, chull=None, diam=None)


def run_one(state, level, group, target_mm, summary_rows, solution_rows):
    print("\n=== %s %s (%s)  target #MM >= %d ===" % (state, level, group, target_mm))

    # ---- build the instance (capture construction failures) ----
    try:
        inst = problem_instance(state, level, False, group)
    except Exception as e:
        print("  CONSTRUCT ERROR:", e)
        traceback.print_exc()
        summary_rows.append(_err_row(state, level, None, None,
                                     "construct error: %s" % e))
        return
    if inst.error_message:
        print("  skipped:", inst.error_message)
        summary_rows.append(_err_row(state, level, None, None,
                                     "skipped: %s" % inst.error_message))
        return

    inst.model_time_limit = PER_MODEL_TIME_LIMIT
    G = inst.G
    n = G.number_of_nodes()
    k = inst.k
    pop, vap, mvap = inst.population, inst.voting_age_population, inst.minority_population
    L, U = inst.L, inst.U
    csr = fg.CSRGraph(n, G.edges())

    # geometry for PP / Reock / convex-hull (None if geopandas/shapefile missing)
    gdf = None
    try:
        import geopandas
        plural = "counties" if level == "county" else "tracts"
        gdf = geopandas.read_file("./raw_data/%s/shape/%s_%s.shp" % (level, state, plural))
    except Exception as e:
        print("  (no geometry scores: %s)" % e)

    s_lo = inst.lower_bound_s
    s_hi = getattr(inst, "upper_bound_s", None)
    try:
        s_hi = int(s_hi)
    except Exception:
        s_hi = n  # fall back to a loose cap

    # The model cannot create more MM districts than the Section-6 upper bound, so
    # cap the target (k_minority==0 => target 0 => the first feasible plan, no
    # wasted climb on instances that admit no MM district at all).
    k_min = getattr(inst, "k_minority", target_mm)
    target_mm = min(target_mm, k_min)

    # Warm-start strategy:
    #  - LARGE instances (tracts, n>300): the exact climb times out with no
    #    incumbent at small s, so solve ONCE at the heuristic plan's diameter d (a
    #    valid s* upper bound) seeded with it -- guarantees a feasible MM incumbent.
    #  - SMALL instances (counties): cold solves are fast, so climb from the lower
    #    bound and CERTIFY the true minimum s* (pinning to the heuristic diameter
    #    would report a suboptimal s, e.g. s=7 when s*=6).
    warm = None
    if n > 300:
        warm = load_gerrychain_warmstart(state, level, G, mvap, vap)
        if warm is not None:
            ws_plan, ws_diam = warm
            s_lo = s_hi = ws_diam
    # small instances (counties) cold-climb from the lower bound for the true s*

    s = s_lo
    total_time = 0.0
    final_s = s
    success = False
    reason = ""

    while True:
        inst.s = s
        final_s = s
        t0 = time.time()
        try:
            inst.labeling_model(True, True, True, bare=False,
                                warm_start=(ws_plan if (warm is not None and s == ws_diam) else None))
        except Exception as e:
            dt = time.time() - t0
            total_time += dt
            print("  s=%d -> EXCEPTION: %s" % (s, e))
            traceback.print_exc()
            reason = "exception at s=%d: %s" % (s, e)
            break

        dt = time.time() - t0
        total_time += dt
        m = inst.m
        nmd = inst.num_minority_districts
        solc = getattr(m, "SolCount", 0)
        status = getattr(m, "status", None)
        objb = _objbound(m)
        # best integer #MM we can read at this s
        if isinstance(nmd, int):
            mm_now = nmd
        elif solc > 0:
            try:
                mm_now = int(round(m.ObjVal))
            except Exception:
                mm_now = None
        else:
            mm_now = None

        print("  s=%d -> num_minority_districts=%s  status=%s  SolCount=%s  "
              "ObjBound=%s  mm_now=%s  (%.1fs, cumulative %.1fs)"
              % (s, nmd, status, solc, objb, mm_now, dt, total_time))

        # ---- decide whether to stop or climb ----
        if isinstance(nmd, int):
            if nmd >= target_mm:
                success = True
                reason = "s*=%d (optimal #MM=%d)" % (s, nmd)
                break
            print("    optimal #MM=%d < target %d at s=%d -> increasing s"
                  % (nmd, target_mm, s))
        elif nmd == "Infeasible":
            print("    infeasible at s=%d -> increasing s" % s)
        else:  # 'Timeout' / 'Model Status: N'
            if solc > 0 and mm_now is not None and mm_now >= target_mm:
                success = True
                reason = "s=%d (timed out, incumbent #MM=%d)" % (s, mm_now)
                break
            if solc > 0:
                print("    %s at s=%d: incumbent #MM=%s < target -> increasing s"
                      % (nmd, s, mm_now))
            else:
                print("    %s at s=%d: NO incumbent found in %ds -> increasing s"
                      % (nmd, s, PER_MODEL_TIME_LIMIT))

        # ---- guards against runaway ----
        if s >= s_hi:
            reason = "hit s upper bound (%d) without reaching #MM>=%d" % (s_hi, target_mm)
            break
        if total_time >= INSTANCE_TIME_BUDGET:
            reason = ("exhausted %.0fs time budget at s=%d without reaching #MM>=%d"
                      % (total_time, s, target_mm))
            break
        s += 1

    print("  STOP: %s" % reason)

    # ---- capture every feasible solution from the final model we solved ----
    try:
        sols = H.extract_all_solutions(inst.m, inst.m._X, inst.m._Y, n, k, csr,
                                       final_s, pop, L, U, mvap, vap, G, "maps",
                                       "%s_%s_ours" % (state, level), gdf=gdf)
    except Exception as e:
        print("  EXTRACT ERROR:", e)
        traceback.print_exc()
        sols = []
    print("  feasible solutions captured:", len(sols))

    # Fallback: if the solve captured nothing but we have a VERIFIED-feasible
    # warm-start plan (it satisfies partition, contiguity, diameter<=s, [L,U] and
    # the MM threshold -- checked by check_feasible_plan), report that plan as our
    # uncertified incumbent rather than leaving the row blank. (Gurobi sometimes
    # discards a borderline start, e.g. an MM district right at 50%.)
    warm_fallback = False
    if not sols and warm is not None:
        ws_minority, ws_majority = ws_plan
        districts = list(ws_minority) + list(ws_majority)
        labels = [None] * n
        for j, D in enumerate(districts):
            for v in D:
                labels[v] = j
        ok, _why = fg.check_feasible_plan(csr, labels, len(districts), final_s,
                                          pop, L, U)
        if ok:
            res = H.score_plan(G, districts, mvap, vap, gdf=gdf, labels=labels)
            H.save_map_csv(G, labels, "maps/%s_%s_ours_warm.csv" % (state, level))
            sols = [dict(sol_index=1, mm=res["mm"], pp=res["avg_pp"],
                         reock=res["avg_reock"], chull=res["avg_chull"],
                         diam=res["diameter"],
                         map_file="maps/%s_%s_ours_warm.csv" % (state, level))]
            warm_fallback = True
            print("  fallback: reporting the verified warm-start plan "
                  "(#MM=%d, diam=%d)" % (res["mm"], res["diameter"]))

    gap = _gap(inst.m)
    best = max(sols, key=lambda r: r["mm"]) if sols else None
    if warm_fallback:
        compact = "s=%d (warm-start plan, not certified)" % final_s
    elif warm is not None:
        compact = "s=%d (GerryChain warm start)" % final_s
    elif success:
        compact = "s*=%d" % final_s
    else:
        compact = "s=%d  [%s]" % (final_s, reason)

    summary_rows.append(dict(
        state=state, level=level, n=n, k=k,
        method="diameter (ours, full)", compact=compact,
        mm=(best["mm"] if best else None),
        time=round(total_time, 1), gap=gap,
        pp=(best["pp"] if best else None),
        reock=(best["reock"] if best else None),
        chull=(best["chull"] if best else None),
        diam=(best["diam"] if best else None),
    ))
    for r in sols:
        solution_rows.append(dict(state=state, level=level,
                                  method="diameter (ours, full)",
                                  gap=gap, time_found=round(total_time, 1), **r))


def _read_existing_rows(path="A1_benchmark_results.xlsx"):
    """Recover (summary_rows, solution_rows) from a previous workbook so a partial
    re-run can MERGE instead of clobber. Returns ([], []) if unreadable."""
    if not os.path.exists(path):
        return [], []
    try:
        from openpyxl import load_workbook
        from benchmark_to_excel import COLS, SOL_COLS
        wb = load_workbook(path, data_only=True)
    except Exception as e:
        print("  (could not read existing workbook: %s)" % e)
        return [], []
    srows, solrows = [], []
    if "A1 benchmark" in wb.sheetnames:
        keys = [k for (k, _h, _f) in COLS]
        for r in wb["A1 benchmark"].iter_rows(min_row=4, values_only=True):
            if any(c is not None for c in r):
                srows.append({keys[j]: r[j] for j in range(min(len(keys), len(r)))})
    if "All feasible solutions" in wb.sheetnames:
        keys = [k for (k, _h, _f) in SOL_COLS]
        for r in wb["All feasible solutions"].iter_rows(min_row=3, values_only=True):
            if any(c is not None for c in r):
                solrows.append({keys[j]: r[j] for j in range(min(len(keys), len(r)))})
    return srows, solrows


if __name__ == "__main__":
    import sys
    wanted = [a.replace("_", ":").upper() for a in sys.argv[1:]]
    instances = INSTANCES
    if wanted:
        instances = [t for t in INSTANCES if ("%s:%s" % (t[0], t[1])).upper() in wanted]
        if not instances:
            print("No instances matched %s. Use e.g.:  MS:county SC:tract"
                  % " ".join(sys.argv[1:]))
            sys.exit(1)

    # when re-running a subset, start from the existing workbook and replace only
    # the rows for the instances we re-run (others are preserved)
    ran = {(t[0], t[1]) for t in instances}
    old_s, old_sol = _read_existing_rows() if wanted else ([], [])
    summary_rows = [r for r in old_s if (r.get("state"), r.get("level")) not in ran]
    solution_rows = [r for r in old_sol if (r.get("state"), r.get("level")) not in ran]

    print("Running %d instance(s): %s"
          % (len(instances), ", ".join("%s %s" % (t[0], t[1]) for t in instances)))
    for (state, level, group, target_mm) in instances:
        try:
            run_one(state, level, group, target_mm, summary_rows, solution_rows)
        except Exception as e:
            print("  FATAL on %s %s: %s" % (state, level, e))
            traceback.print_exc()
            summary_rows.append(_err_row(state, level, None, None,
                                         "fatal: %s" % e))
    # keep instance order stable in the workbook
    order = {"%s:%s" % (t[0], t[1]): i for i, t in enumerate(INSTANCES)}
    summary_rows.sort(key=lambda r: order.get("%s:%s" % (r.get("state"), r.get("level")), 999))
    build_workbook(summary_rows, solution_rows, "A1_benchmark_results.xlsx")
    print("\nWrote A1_benchmark_results.xlsx")
