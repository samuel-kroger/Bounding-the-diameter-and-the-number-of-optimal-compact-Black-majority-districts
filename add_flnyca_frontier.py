"""
add_flnyca_frontier.py -- add the FL, NY and CA rows to the Figure 1 data.

WHY NO RE-RUN IS NEEDED
-----------------------
Figure 1 is drawn from results_tradeoff/tradeoff_frontier.csv, whose
"frontier s_min(k_m)" column lists, for each representation level k_m, the
smallest diameter at which a plan with that many majority-minority districts was
found. Those values are heuristic incumbents, i.e. UPPER bounds on s*(k_m), which
is exactly what the GerryChain sweep already produced and recorded in
gerrychain_FLNYCA.csv as its Pareto frontier.

So the new rows can be derived directly, with one observation: a plan with m
majority-minority districts is simultaneously a plan for every level below m.
Monotonicity therefore fills in the intermediate levels. For California the
sweep returned (20|11) and (21|14), which yields

    k_m =  1..11  ->  s = 20
    k_m = 12..14  ->  s = 21

This is the same reasoning the existing rows use, so the new rows are consistent
with the ones already plotted.

USAGE
    python add_flnyca_frontier.py            # preview, writes nothing
    python add_flnyca_frontier.py --write    # append the rows
"""
import argparse
import csv
import os
import re
import shutil

FRONTIER = "results_tradeoff/tradeoff_frontier.csv"
GC = "gerrychain_FLNYCA.csv"
NEW_STATES = ("FL", "NY", "CA")

# k, from the apportionment table
K = {"FL": 28, "NY": 26, "CA": 52}

# k_m upper bounds produced by the sweep (Table 3 column "kmUB(Tab3)")
KM_UB = {("FL", "black"): 4,  ("FL", "hispanic"): 10,
         ("NY", "black"): 6,  ("NY", "hispanic"): 5,
         ("CA", "black"): 0,  ("CA", "hispanic"): 33}


def parse_frontier(text):
    """'(21|14) (20|11)' -> {14: 21, 11: 20}"""
    out = {}
    for s, mm in re.findall(r"\((\d+)\|(\d+)\)", text or ""):
        s, mm = int(s), int(mm)
        if mm and (mm not in out or s < out[mm]):
            out[mm] = s
    return out


def expand(points):
    """Fill intermediate levels by monotonicity: a plan with m districts also
    serves every level below m, at the same diameter."""
    if not points:
        return {}, 0
    best = max(points)
    filled = {}
    for level in range(1, best + 1):
        # smallest diameter among plans achieving at least this level
        cands = [s for mm, s in points.items() if mm >= level]
        filled[level] = min(cands)
    return filled, best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(GC):
        raise SystemExit("missing %s -- run merge_gc.py first" % GC)

    with open(GC) as fh:
        gc_rows = list(csv.DictReader(fh))

    new_rows = []
    for r in gc_rows:
        st, grp = r["State"].strip(), r["group"].strip()
        if st not in NEW_STATES:
            continue
        pts = parse_frontier(r.get("pareto frontier", ""))
        filled, best = expand(pts)
        frontier_str = " ".join("%d:%d" % (lvl, filled[lvl])
                                for lvl in sorted(filled))
        new_rows.append([st, "tract", grp, K[st], KM_UB[(st, grp)],
                         best, frontier_str, ""])

    order = {"black": 0, "hispanic": 1}
    new_rows.sort(key=lambda x: (order[x[2]], x[0]))

    print("Rows to add:\n")
    hdr = ["State", "level", "group", "k", "kmUB(Tab3)", "max_MM_found",
           "frontier s_min(k_m)", "s*_proven"]
    print("  " + " | ".join(hdr))
    for r in new_rows:
        print("  " + " | ".join(str(c) for c in r))

    if not args.write:
        print("\n(preview only; pass --write to append)")
        return

    with open(FRONTIER) as fh:
        existing = [l.rstrip("\n") for l in fh if l.strip()]
    kept = [existing[0]] + [l for l in existing[1:]
                            if l.split(",")[0] not in NEW_STATES]
    shutil.copy(FRONTIER, FRONTIER + ".bak")
    with open(FRONTIER, "w", newline="") as fh:
        fh.write("\n".join(kept) + "\n")
        w = csv.writer(fh)
        w.writerows(new_rows)
    print("\nWrote %s (backup at %s.bak)" % (FRONTIER, FRONTIER))


if __name__ == "__main__":
    main()
