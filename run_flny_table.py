"""
run_flny_table.py -- run one paper table for the FL, NY and CA tract instances
WITHOUT clobbering the results already on disk for the paper's instances.

WHY THIS WRAPPER EXISTS
-----------------------
run_table.py deliberately deletes a table's CSV before running, so that a re-run
produces a clean table instead of appending duplicate rows. That is the right
behavior when you are regenerating the whole paper, but it is destructive here:
running the FL/NY sweep directly would wipe results_black/instance_info.csv and
friends, i.e. the numbers behind the tables already in the manuscript.

So this wrapper:
  1. snapshots the target CSV(s),
  2. runs run_table.py for the FL/NY dataset,
  3. harvests the FL/NY output to <name>_FL_NY.csv,
  4. restores the original CSV(s) byte for byte.

The per-state caches that comp_experiment writes (results_<group>/lower_bound_s/
<STATE>_<level>.txt and upper_bound_minority_districts/...) are keyed by state,
so FL and NY only ADD files there. Nothing existing is touched.

USAGE
    python run_flny_table.py --table lower_bound_s --group hispanic
    python run_flny_table.py --table lower_bound_s --group black
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

# dataset key in data.json for each group
DATASET = {"black": "flnyca", "hispanic": "flnyca_hispanic"}

# CSVs each table writes (mirrors run_table.TABLE_CSV, plus the resumable ones
# that run_table intentionally does not wipe but which we still must protect,
# because the FL/NY rows would otherwise be appended into the paper's file).
# Paths starting with "@" are repo-root files rather than results_<group> files;
# generate_gerrychain_run writes gerrychain.csv at the root.
TABLE_CSV = {
    "instance_info":     ["instance_info.csv"],
    "lower_bound_s":     ["lower_bound_s_table.csv"],
    "upper_bound_mm":    ["upper_bound_minority_table.csv"],
    "fixing":            ["fixing_info.csv", "fixing_info_each_pass.csv"],
    "symmetry":          ["symmetry_table.csv"],
    "symmetry_feasible": ["symmetry_feasible_table.csv"],
    "gerrychain":        ["@gerrychain.csv"],
    "point_check":       ["point_check.csv"],
    "upper_bound_s":     ["upper_bound_s_table.csv"],
    "interesting_s":     ["interesting_s_vals.csv"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, choices=sorted(TABLE_CSV))
    ap.add_argument("--group", required=True, choices=["black", "hispanic"])
    args = ap.parse_args()

    resdir = "results_%s" % args.group
    os.makedirs(resdir, exist_ok=True)

    # Pre-flight: problem_instance writes its per-state caches into these
    # subdirectories but never creates them, so a state that has never been run
    # before (FL, NY) dies with FileNotFoundError partway through __init__.
    # Creating them up front is harmless for states already on disk.
    for sub in ("fixing", "fixing/txt", "gerrychain", "interesting_s_vals",
                "lower_bound_s", "max_independent_set",
                "upper_bound_minority_districts", "upper_bound_s"):
        os.makedirs(os.path.join(resdir, sub), exist_ok=True)

    targets = [fn[1:] if fn.startswith("@") else os.path.join(resdir, fn)
               for fn in TABLE_CSV[args.table]]

    # 1. snapshot
    stash = tempfile.mkdtemp(prefix="flny_stash_")
    saved = {}
    for t in targets:
        if os.path.exists(t):
            dst = os.path.join(stash, os.path.basename(t))
            shutil.copy2(t, dst)
            saved[t] = dst
            print("[flny] protected %s" % t)
        # remove so the run starts clean even for the resumable tables
        if os.path.exists(t):
            os.remove(t)

    # 2. run
    cmd = [sys.executable, "run_table.py",
           "--table", args.table,
           "--group", args.group,
           "--dataset", DATASET[args.group]]
    print("[flny] %s" % " ".join(cmd))
    sys.stdout.flush()
    rc = subprocess.call(cmd)

    # 3. harvest, then 4. restore
    try:
        for t in targets:
            if os.path.exists(t):
                # include the group: root-level outputs (gerrychain.csv) are not
                # separated by a results_<group> directory, so without this the
                # second group's harvest silently overwrites the first's
                suffix = "_FLNYCA_%s.csv" % args.group
                out = t.replace(".csv", suffix)
                shutil.move(t, out)
                print("[flny] FL/NY results -> %s" % out)
            else:
                print("[flny] NOTE: %s produced no rows" % t)
    finally:
        for t, src in saved.items():
            shutil.copy2(src, t)
            print("[flny] restored %s" % t)
        shutil.rmtree(stash, ignore_errors=True)

    sys.exit(rc)


if __name__ == "__main__":
    main()
