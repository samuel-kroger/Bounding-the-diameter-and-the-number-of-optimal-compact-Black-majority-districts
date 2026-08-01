"""
score_belotti_plans.py -- re-score Belotti/Buchanan/Ezazipour's exported plans
with OUR scorer, so PP / Reock / convex-hull / diameter / #MM come from the same
geometry code we use on our own plans (apples-to-apples for the reviewer).

It walks their results folder for the GEOID20,district plan CSVs written by
run_belotti_one.py (named <state>_<level>_<objective>_<contiguity>.csv), maps
each plan onto OUR graph + shapefile by GEOID20, runs
benchmark_run_helpers.score_plan, and writes the numbers into the workbook's
"Belotti et al. (published, PP-optimal)" rows -- merged side-by-side with our
"diameter (ours, full)" rows from A1_benchmark_results.xlsx when that file is
present.

Run from the code/ directory AFTER run_belotti_our_instances.bat finishes:
    python score_belotti_plans.py [path\to\results_<date>]
If no path is given it auto-finds the newest ../Polsby_Popper_optimization-main/results_*.
"""
import csv
import glob
import json
import os
import sys

import benchmark_run_helpers as H
from benchmark_to_excel import build_workbook, COLS

try:
    import networkx as nx
    from networkx.readwrite import json_graph
except Exception as e:  # pragma: no cover
    print("networkx is required:", e)
    sys.exit(1)

PP_METHOD = "Belotti et al. (PP-optimal)"
CCC_METHOD = "Belotti et al. (CCC, MM-seeking)"
GC_METHOD = "GerryChain (heuristic UB)"
OUR_METHOD = "diameter (ours, full)"


def method_for(objective):
    if objective == "gingles":
        return CCC_METHOD
    if objective == "gerrychain":
        return GC_METHOD
    return PP_METHOD


# our instance set -> minority group (selects the VAP race field)
GROUP = {
    ("MS", "county"): "black", ("AL", "county"): "black",
    ("LA", "county"): "black", ("SC", "county"): "black",
    ("NM", "tract"): "hispanic", ("MS", "tract"): "black",
    ("SC", "tract"): "black", ("LA", "tract"): "black",
    ("CO", "tract"): "hispanic",
}
MIN_FIELD = {"black": "P0030004", "hispanic": "P0040002"}  # "alone" (our model's def)
VAP_FIELD = "P0030001"
# "any-part" Black VAP codes (alone OR in combination) -- their CCC carve uses
# this definition. Hispanic any-part == P0040002, same as the alone field.
BLACK_ANYPART = [
    "P0030004", "P0030011", "P0030016", "P0030017", "P0030018", "P0030019",
    "P0030027", "P0030028", "P0030029", "P0030030", "P0030037", "P0030038",
    "P0030039", "P0030040", "P0030041", "P0030042", "P0030048", "P0030049",
    "P0030050", "P0030051", "P0030052", "P0030053", "P0030058", "P0030059",
    "P0030060", "P0030061", "P0030064", "P0030065", "P0030066", "P0030067",
    "P0030069", "P0030071",
]


def find_results_dirs():
    """ALL results_* folders (so invpp, gingles, etc. -- possibly written on
    different days -- are all picked up)."""
    if len(sys.argv) > 1:
        return [sys.argv[1]]
    cands = glob.glob(os.path.join("..", "Polsby_Popper_optimization-main", "results_*"))
    cands = [c for c in cands if os.path.isdir(c)]
    if not cands:
        print("No results_* folder found. Pass the path explicitly:")
        print("    python score_belotti_plans.py <path-to-results_DATE>")
        sys.exit(1)
    return sorted(cands, key=os.path.getmtime)


def load_instance(state, level):
    """Our graph + per-node vap/minority/geoid (indexed by node id), and a gdf."""
    plural = "counties" if level == "county" else "tracts"
    with open("raw_data/%s/json/%s_%s.json" % (level, state, plural)) as f:
        data = json.load(f)
    G = json_graph.adjacency_graph(data, multigraph=False)
    n = G.number_of_nodes()
    group = GROUP[(state, level)]
    minkey = MIN_FIELD[group]
    vap = [0] * n
    mino = [0] * n           # "alone" def (matches our diameter model)
    mino_any = [0] * n       # "any-part" def (matches their CCC carve)
    geoid_to_idx = {}
    for nd in data["nodes"]:
        i = nd["id"]
        vap[i] = nd[VAP_FIELD]
        mino[i] = nd[minkey]
        if group == "black":
            mino_any[i] = sum(nd[c] for c in BLACK_ANYPART)
        else:  # hispanic any-part == alone field
            mino_any[i] = nd[minkey]
        g = str(nd["GEOID20"]).strip()
        geoid_to_idx[g] = i
        geoid_to_idx[g.zfill(5 if level == "county" else 11)] = i
    gdf = None
    try:
        import geopandas
        gdf = geopandas.read_file("raw_data/%s/shape/%s_%s.shp" % (level, state, plural))
    except Exception as e:
        print("    (no geometry for %s %s: %s)" % (state, level, e))
    return G, vap, mino, mino_any, geoid_to_idx, gdf


def labels_from_their_csv(path, geoid_to_idx, n, level):
    labels = [None] * n
    missing = 0
    for r in csv.DictReader(open(path)):
        g = str(r["GEOID20"]).strip()
        idx = geoid_to_idx.get(g) or geoid_to_idx.get(g.zfill(5 if level == "county" else 11))
        if idx is None:
            missing += 1
            continue
        labels[idx] = int(r["district"])
    return labels, missing


def read_their_summary(results_dirs):
    """(state,level,objective) -> (MIP_time, gap) from every summary.csv."""
    out = {}
    for results_dir in results_dirs:
        p = os.path.join(results_dir, "summary.csv")
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p)):
            try:
                t = float(r["MIP_time"])
            except Exception:
                t = None
            gap = None
            try:
                obj = float(r["MIP_obj"]); bnd = float(r["MIP_bound"])
                if obj:
                    gap = abs(obj - bnd) / abs(obj)
            except Exception:
                pass
            out[(r["state"], r["level"], r["objective"])] = (t, gap)
    return out


def score_their_plans(results_dirs):
    meta = read_their_summary(results_dirs)
    summary_by_key, solution_rows = {}, []
    plan_csvs = []
    for d in results_dirs:
        plan_csvs += sorted(glob.glob(os.path.join(d, "*_*_*_*.csv")))
    for path in plan_csvs:
        base = os.path.splitext(os.path.basename(path))[0]
        parts = base.split("_")
        if len(parts) < 4 or base == "summary":
            continue
        state, level, objective, contiguity = parts[0], parts[1], parts[2], parts[3]
        if (state, level) not in GROUP:
            print("  [skip] %s (not in our instance set)" % base)
            continue
        method = method_for(objective)
        print("  scoring %s  ->  %s" % (base, method))
        G, vap, mino, mino_any, geoid_to_idx, gdf = load_instance(state, level)
        n = G.number_of_nodes()
        labels, missing = labels_from_their_csv(path, geoid_to_idx, n, level)
        if missing:
            print("    WARNING: %d GEOIDs in their plan did not match our graph; "
                  "skipping geometry for safety" % missing)
            gdf = None
        if any(l is None for l in labels):
            gdf = None
        k = max((l for l in labels if l is not None), default=-1) + 1
        districts = H.districts_from_labels(labels, k)
        # headline #MM uses the ANY-PART definition (VRA standard, matching both
        # our re-defined model and their CCC carve). Black-alone is reported in
        # the note for transparency.
        res = H.score_plan(G, districts, mino_any, vap, gdf=gdf, labels=labels)
        mm_any = res["mm"]
        mm_alone = H.num_minority_districts(G, districts, mino, vap)
        t, gap = meta.get((state, level, objective), (None, None))
        if objective == "gingles":
            tag = "CCC MVAP>=0.5"
        elif objective == "gerrychain":
            tag = "GerryChain UB"
        else:
            tag = "PP-opt(%s)" % objective
        compact = "%s | MMalone=%d" % (tag, mm_alone)
        # dedup by (state, level, method): keep the most recently scored
        summary_by_key[(state, level, method)] = dict(
            state=state, level=level, n=n, k=k, method=method, compact=compact,
            mm=mm_any, time=t, gap=gap,
            pp=res["avg_pp"], reock=res["avg_reock"], chull=res["avg_chull"],
            diam=res["diameter"])
        solution_rows.append(dict(
            state=state, level=level, method=method, sol_index=1,
            mm=mm_any, time_found=t, gap=gap,
            pp=res["avg_pp"], reock=res["avg_reock"], chull=res["avg_chull"],
            diam=res["diameter"], map_file=os.path.relpath(path)))
    return list(summary_by_key.values()), solution_rows


def read_our_summary(path="A1_benchmark_results.xlsx"):
    """Best-effort recovery of our 'diameter (ours, full)' rows from the workbook
    so we can merge both methods. Returns [] if the file is unreadable."""
    if not os.path.exists(path):
        return []
    try:
        from openpyxl import load_workbook
        wb = load_workbook(path, data_only=True)
        ws = wb["A1 benchmark"]
    except Exception as e:
        print("  (could not read our rows from %s: %s)" % (path, e))
        return []
    keys = [k for (k, _h, _f) in COLS]
    rows = []
    for r in ws.iter_rows(min_row=4, values_only=True):
        if not any(c is not None for c in r):
            continue
        row = {keys[j]: r[j] for j in range(min(len(keys), len(r)))}
        if row.get("method") == OUR_METHOD:
            rows.append(row)
    return rows


def merge(our_rows, their_rows):
    """Group by (state, level): our row first, then their method rows
    (PP-optimal, then CCC), preserving the order instances first appear."""
    order = []
    seen = set()
    by_key = {}
    for row in our_rows + their_rows:
        key = (row.get("state"), row.get("level"))
        if key not in seen:
            seen.add(key)
            order.append(key)
        by_key.setdefault(key, {"ours": None, "theirs": []})
    for row in our_rows:
        by_key[(row["state"], row["level"])]["ours"] = row
    for row in their_rows:
        by_key[(row["state"], row["level"])]["theirs"].append(row)
    rank = {OUR_METHOD: 0, GC_METHOD: 1, PP_METHOD: 2, CCC_METHOD: 3}
    merged = []
    for key in order:
        grp = by_key[key]
        if grp["ours"]:
            merged.append(grp["ours"])
        for row in sorted(grp["theirs"], key=lambda r: rank.get(r["method"], 9)):
            merged.append(row)
    return merged


def main():
    results_dirs = find_results_dirs()
    print("Scoring Belotti plans in:", ", ".join(results_dirs))
    their_rows, their_sols = score_their_plans(results_dirs)
    if not their_rows:
        print("No matching plan CSVs found. Did the batch run finish?")
        return
    our_rows = read_our_summary()
    merged = merge(our_rows, their_rows)
    build_workbook(merged, their_sols, "A1_benchmark_combined.xlsx")
    print("\nWrote A1_benchmark_combined.xlsx  (%d ours + %d theirs rows)."
          % (len(our_rows), len(their_rows)))
    # console preview ("#MM" is the alone-definition count; any-part is in note)
    print("\n%-6s %-7s %-34s %4s %7s %7s %7s %5s  %s"
          % ("state", "level", "method", "#MM", "PP", "Reock", "CH", "diam", "note"))
    for row in merged:
        print("%-6s %-7s %-34s %4s %7s %7s %7s %5s  %s" % (
            row.get("state"), row.get("level"), row.get("method"),
            row.get("mm"),
            ("%.3f" % row["pp"]) if row.get("pp") is not None else "-",
            ("%.3f" % row["reock"]) if row.get("reock") is not None else "-",
            ("%.3f" % row["chull"]) if row.get("chull") is not None else "-",
            row.get("diam"),
            row.get("compact") or ""))


if __name__ == "__main__":
    main()
