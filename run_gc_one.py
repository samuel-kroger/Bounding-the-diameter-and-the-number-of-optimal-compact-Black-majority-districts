"""
run_gc_one.py -- run the GerryChain heuristic for ONE (state, group) pair and
write its own CSV, so several instances can run in parallel.

WHY
---
GerryChain is single threaded and is by far the longest step in the sweep
(roughly 2.3 s per iteration on a tract instance, so ~6.4 h for the standard
10,000 iterations). Run sequentially, six instances take about 38 hours. Run in
parallel on a machine with plenty of cores they take about 6.4 hours, because
they do not compete for anything.

The paper's own driver writes every instance into a single shared gerrychain.csv,
which cannot be done safely from parallel processes. This script gives each
instance a private output file; merge_gc.py stitches them back together.

USAGE
    python run_gc_one.py --state FL --group hispanic
    python run_gc_one.py --state CA --group black --iters 10000

OUTPUT
    gerrychain_<STATE>_<group>.csv        (skipped if it already exists)
    gc_log_<STATE>_<group>.txt
"""
import argparse
import os
import sys
import time
import traceback


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--group", required=True, choices=["black", "hispanic"])
    ap.add_argument("--level", default="tract")
    ap.add_argument("--iters", type=int,
                    default=int(os.environ.get("GC_ITERS", 10000)))
    ap.add_argument("--force", action="store_true",
                    help="re-run even if the output already exists")
    args = ap.parse_args()

    tag = "%s_%s" % (args.state, args.group)
    out = "gerrychain_%s.csv" % tag

    # resumability: an instance that already produced output is not redone, so
    # an interrupted sweep can simply be restarted
    if os.path.exists(out) and not args.force:
        print("[gc] %s already done (%s); skipping" % (tag, out))
        return 0

    sys.path.insert(0, os.getcwd())
    from comp_experiment import problem_instance

    t0 = time.time()
    print("[gc] %s: %d iterations, starting %s"
          % (tag, args.iters, time.strftime("%H:%M:%S")), flush=True)

    try:
        inst = problem_instance(args.state, args.level, False, args.group)
        if getattr(inst, "error_message", ""):
            print("[gc] %s SKIPPED: %s" % (tag, inst.error_message))
            return 0

        results = inst.run_GerryChain_heuristic(args.iters)
        elapsed = time.time() - t0

        # Persist each frontier plan as a .pckl, exactly as
        # generate_gerrychain_run does. This is ESSENTIAL and was missing in the
        # first version of this script: check_points_with_prescribed_model calls
        # _load_heuristic_plans(), which globs
        #     results_<group>/gerrychain/<STATE>_<level>_(ss,mm).pckl
        # and uses those plans to CERTIFY feasibility without solving a MIP.
        # Writing only the summary CSV left the frontier step with nothing to
        # warm start from, which is why it timed out on FL and CA even after the
        # chains had found good plans.
        import pickle
        outdir = './results_%s/gerrychain/' % args.group
        os.makedirs(outdir, exist_ok=True)
        saved = 0
        for result in results:
            partition, gc_s, gc_mm = result[0], result[1], result[2]
            minority_districts, majority_districts = [], []
            for district in range(inst.k):
                nodes, pop, minority = [], 0, 0
                for node in partition.graph:
                    if partition.assignment[node] == district:
                        nodes.append(node)
                        pop += inst.voting_age_population[node]
                        minority += inst.minority_population[node]
                if minority > inst.f * pop:
                    minority_districts.append(nodes)
                else:
                    majority_districts.append(nodes)
            fn = '%s%s_%s_(%02d,%02d).pckl' % (outdir, args.state, args.level,
                                               gc_s, gc_mm)
            with open(fn, 'wb') as fh:
                pickle.dump([minority_districts, majority_districts], fh)
            saved += 1
        print("[gc] %s saved %d plan file(s) for the frontier step"
              % (tag, saved), flush=True)

        # run_GerryChain_heuristic returns the Pareto frontier as a list of
        # [partition, s, #MM] triples (NOT (s, #MM) pairs), matching how
        # generate_gerrychain_run reads result[1] and result[2].
        frontier = " ".join("(%s|%s)" % (r[1], r[2]) for r in results) \
            if results else ""
        best = max((r[2] for r in results), default=0) if results else 0
        best_s = min((r[1] for r in results if r[2] == best), default="") \
            if results else ""

        with open(out, "w") as fh:
            fh.write("State,parcel level,group,iterations,time,best #MM,"
                     "s at best,pareto frontier\n")
            fh.write("%s,%s,%s,%d,%.2f,%d,%s,%s\n"
                     % (args.state, args.level, args.group, args.iters,
                        elapsed, best, best_s, frontier))
        print("[gc] %s DONE in %.1f h, best #MM = %d at s = %s, frontier = %s"
              % (tag, elapsed / 3600.0, best, best_s, frontier or "(none)"),
              flush=True)
        return 0

    except Exception as exc:                                   # noqa: BLE001
        print("[gc] %s FAILED after %.1f min: %s"
              % (tag, (time.time() - t0) / 60.0, exc), flush=True)
        traceback.print_exc()
        with open("gc_error_%s.log" % tag, "w") as fh:
            fh.write("%s\n\n" % exc)
            traceback.print_exc(file=fh)
        return 1


if __name__ == "__main__":
    sys.exit(main())
