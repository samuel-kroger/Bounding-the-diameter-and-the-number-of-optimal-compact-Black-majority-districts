"""
run_table.py -- reproducible, parameterized driver for the paper's tables.

This REPLACES the old workflow of hand-editing comp_experiment.py __main__ (setting
`group`/`dataset` and uncommenting one method). Every table is now produced by a
single command, so the whole suite is reproducible on one machine.

Usage:
    python run_table.py --table <name> --group <black|hispanic|none> --dataset <name>

Tables (the comp_experiment method each dispatches to, and the CSV it writes under
results_<group>/):
    instance_info   -> generate_instance_info             (instance_info.csv)
    lower_bound_s   -> generate_lower_bound_s_table        (lower_bound_s_table.csv)
    upper_bound_mm  -> generate_upper_bound_minority_table (upper_bound_minority_table.csv)
    fixing          -> generate_fixing_info (+ aggregate)  (fixing_info.csv)
    symmetry        -> generate_symmetry_table             (symmetry_table.csv)
    point_check     -> check_points_with_prescribed_model  (point_check.csv)   [s vs #MM]
    upper_bound_s   -> create_upper_bound_s_table          (upper_bound_s_table.csv)
    interesting_s   -> find_interesting_s_vals_no_minority (interesting_s_vals.csv)

The heavy numerical paths these methods rely on (the length-s separator callback,
the three-phase fixing pre-screen, the independence-number / power-graph routines,
and the up-front distance cuts) all use the array-based fast_graph routines rather
than networkx; see fast_graph.py.

Examples:
    python run_table.py --table lower_bound_s  --group black    --dataset paper
    python run_table.py --table lower_bound_s  --group hispanic --dataset hispanic_paper
"""
import argparse
import json
import os
import sys
import time
import traceback

TABLES = {
    "instance_info":  "generate_instance_info",
    "lower_bound_s":  "generate_lower_bound_s_table",
    "upper_bound_mm": "generate_upper_bound_minority_table",
    "fixing":         "generate_fixing_info",
    "symmetry":          "generate_symmetry_table",           # Panel A: s = ell_s (infeasible regime)
    "symmetry_feasible": "generate_symmetry_feasible_table",  # Panel B: s from data.json (feasible regime)
    "gerrychain":     "generate_gerrychain_run",             # heuristic plans (feasibility certificates for Table 6)
    "point_check":    "check_points_with_prescribed_model",
    "upper_bound_s":  "create_upper_bound_s_table",
    "interesting_s":  "find_interesting_s_vals_no_minority",
}

# the CSV(s) each table writes; removed at the start so a re-run produces a CLEAN
# table instead of appending to (and duplicating) a previous run's rows.
TABLE_CSV = {
    "instance_info":  ["instance_info.csv"],
    "lower_bound_s":  ["lower_bound_s_table.csv"],
    "upper_bound_mm": ["upper_bound_minority_table.csv"],
    "fixing":         ["fixing_info.csv", "fixing_info_each_pass.csv"],
    # NB: "symmetry" is intentionally omitted -- generate_symmetry_table is
    # RESUMABLE (it skips (state, level) rows already present), so its CSV must
    # NOT be wiped at the start of a run. That lets an interrupted ~day-long
    # sweep be continued by simply re-running the batch. To start the symmetry
    # table fresh, delete results_<group>/symmetry_table.csv by hand.
    # NB: "point_check" omitted -- the frontier sweep is RESUMABLE (it skips
    # (state, level, k_m) rows already recorded and restarts the monotone carry
    # from the largest feasible s* found), so its CSV must not be wiped.
    # Delete results_<group>/point_check.csv by hand to start it fresh.
    "upper_bound_s":  ["upper_bound_s_table.csv"],
    "interesting_s":  ["interesting_s_vals.csv"],
}


def main():
    ap = argparse.ArgumentParser(description="Regenerate one paper table.")
    ap.add_argument("--table", required=True, choices=sorted(TABLES))
    ap.add_argument("--group", required=True, choices=["black", "hispanic", "none"])
    ap.add_argument("--dataset", required=True, help="a key in data.json (e.g. paper, hispanic_paper, county)")
    args = ap.parse_args()

    with open("data.json") as f:
        data = json.load(f)
    if args.dataset not in data:
        print("Unknown dataset '%s'. Choices: %s" % (args.dataset, ", ".join(sorted(data))))
        sys.exit(1)
    reqs = data[args.dataset]

    # import here so --help works without gurobipy installed
    from comp_experiment import problem_instance

    method = TABLES[args.table]
    os.makedirs("results_%s" % args.group, exist_ok=True)
    # start clean: remove this table's CSV(s) so the run does not append to / dup
    # a previous run's rows
    for fn in TABLE_CSV.get(args.table, []):
        p = os.path.join("results_%s" % args.group, fn)
        if os.path.exists(p):
            os.remove(p)
            print("  (cleared stale %s)" % p)
    print("=== Table '%s' | method %s() | group=%s | dataset=%s | %d instance(s) ==="
          % (args.table, method, args.group, args.dataset, len(reqs)))

    t_all = time.time()
    last_inst = None
    for req in reqs:
        st, lvl = req["state"], req["parcel_level"]
        print("\n--- %s %s ---" % (st, lvl))
        sys.stdout.flush()
        try:
            inst = problem_instance(st, lvl, req.get("s"), args.group)
            if getattr(inst, "error_message", ""):
                print("  skipped:", inst.error_message)
                continue
            getattr(inst, method)()
            last_inst = inst
        except Exception as e:
            print("  ERROR on %s %s: %s" % (st, lvl, e))
            traceback.print_exc()
            # also persist the traceback so a failure during an unattended,
            # console-less overnight run is not lost (this is how we diagnose
            # instances that silently produce no row).
            try:
                logp = os.path.join("results_%s" % args.group, "%s_errors.log" % args.table)
                with open(logp, "a") as lf:
                    lf.write("\n=== ERROR on %s %s at %s ===\n"
                             % (st, lvl, time.strftime("%Y-%m-%d %H:%M:%S")))
                    traceback.print_exc(file=lf)
            except Exception:
                pass

    # the fixing table is written per pass; aggregate it once at the end
    if args.table == "fixing" and last_inst is not None:
        try:
            last_inst.create_fixing_csv_file()
        except Exception as e:
            print("  (fixing aggregation failed: %s)" % e)

    print("\nDone in %.1fs. Output: results_%s/  (see the table's CSV)."
          % (time.time() - t_all, args.group))


if __name__ == "__main__":
    main()
