"""
merge_gc.py -- stitch the per-instance gerrychain_<STATE>_<group>.csv files
written by the parallel step into one table.

USAGE
    python merge_gc.py                    # writes gerrychain_FLNYCA.csv
"""
import csv
import glob
import os

STATES = ("FL", "NY", "CA")
OUT = "gerrychain_FLNYCA.csv"
HEADER = ["State", "parcel level", "group", "iterations", "time",
          "best #MM", "s at best", "pareto frontier"]


def main():
    rows = []
    # Match ONLY the per-instance files this sweep writes, naming the states
    # explicitly. A looser glob picks up leftovers from earlier runs
    # (gerrychain_FL_NY.csv) and any scratch files, which either have a
    # different column layout or are not real instances.
    paths = sorted(set(
        p for st in STATES for g in ("black", "hispanic")
        for p in glob.glob("gerrychain_%s_%s.csv" % (st, g))))
    for path in paths:
        if os.path.basename(path) == OUT:
            continue
        with open(path) as fh:
            r = [row for row in csv.reader(fh) if row]
        if len(r) < 2:
            print("  ! %s has no data row, skipping" % path)
            continue
        for row in r[1:]:
            if len(row) != len(HEADER):
                print("  ! %s: expected %d columns, got %d, skipping"
                      % (path, len(HEADER), len(row)))
                continue
            rows.append(row)
        print("  + %s" % path)

    if not rows:
        print("No per-instance GerryChain files found. Nothing to merge.")
        return

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        w.writerows(rows)

    print("\nWrote %s with %d row(s):\n" % (OUT, len(rows)))
    widths = [max(len(str(x)) for x in [HEADER[i]] + [r[i] for r in rows])
              for i in range(len(HEADER))]
    print("  " + "  ".join(h.ljust(widths[i]) for i, h in enumerate(HEADER)))
    for r in rows:
        print("  " + "  ".join(str(c).ljust(widths[i]) for i, c in enumerate(r)))


if __name__ == "__main__":
    main()
