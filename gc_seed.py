"""
gc_seed.py -- time-boxed, retrying seed generation for GerryChain.

THE PROBLEM
-----------
GerryChain builds its initial plan with recursive_tree_part, which repeatedly
draws a random spanning tree and looks for a population-balanced edge to cut.
At the +/-0.5% deviation this paper uses, that search can fail. GerryChain's
default budget is

    method = partial(bipartition_tree, max_attempts=10000)

and when it runs out it raises

    RuntimeError: Could not find a possible cut after 10000 attempts

killing the instance before the Markov chain takes a single step. This is what
happened to Florida on the first sweep.

Two things are worth being clear about:

  * This is a SEEDING budget, not the chain length. Raising the number of
    GerryChain iterations cannot help, because the failure happens before
    iteration 1. (Raising iterations does help RESULT QUALITY once seeding
    succeeds; that knob is GC_ITERS.)

  * The failure is stochastic AND heavy-tailed. Measured here: FL and NY seed in
    3 to 7 seconds on most random seeds, but some seeds run for many minutes
    without succeeding even with a 200,000-attempt budget. Simply raising
    max_attempts therefore trades a fast crash for a slow hang, which is worse
    for an unattended overnight sweep.

THE FIX
-------
Give each seed attempt a wall-clock budget in a separate process. If it does not
finish in time, kill it and retry from a different random seed. Because the
runtime distribution is heavy-tailed, restarting is far more effective than
waiting: several short attempts beat one long one.

The work runs in a subprocess because recursive_tree_part spends its time in
tight loops that ignore thread-level interruption, so a thread could not be
reclaimed. A spawned process can simply be terminated.

TUNING (environment variables)
    GC_SEED_TIMEOUT    seconds per attempt          (default 120)
    GC_SEED_RETRIES    number of attempts           (default 40)
    GC_SEED_ATTEMPTS   bipartition_tree max_attempts(default 25000)
    GC_SEED_RESELECT   1 to allow pair reselection  (default 0)

Defaults give up to 40 x 120 s = 80 minutes of seeding effort per instance,
which in testing was far more than any of FL, NY or CA needed.
"""
import multiprocessing as mp
import os
import random
import time
import warnings


def _worker(json_path, k, pop_col, epsilon, max_attempts, reselect, seed, q):
    """Runs in a child process; puts the assignment dict on the queue."""
    warnings.filterwarnings("ignore")
    try:
        from functools import partial
        from gerrychain import Graph
        from gerrychain.tree import recursive_tree_part, bipartition_tree

        G = Graph.from_json(json_path)
        target = sum(G.nodes[i][pop_col] for i in G.nodes()) / k
        method = partial(bipartition_tree,
                         max_attempts=max_attempts,
                         allow_pair_reselection=reselect)
        random.seed(seed)
        assignment = recursive_tree_part(G, range(k), target, pop_col,
                                         epsilon, 1, method=method)
        # keys may be numpy ints; normalize so the parent can use them directly
        q.put(("ok", {int(v): int(d) for v, d in assignment.items()}))
    except Exception as exc:                                   # noqa: BLE001
        q.put(("err", "%s: %s" % (type(exc).__name__, exc)))


def seed_partition(json_path, k, pop_col="P0010001", epsilon=0.005,
                   timeout=None, retries=None, max_attempts=None,
                   reselect=None, log=print):
    """Return an assignment dict, or None if every attempt failed.

    Each attempt runs in its own process with a wall-clock timeout; on timeout
    or error we retry from a fresh random seed.
    """
    timeout = int(os.environ.get("GC_SEED_TIMEOUT", timeout or 120))
    retries = int(os.environ.get("GC_SEED_RETRIES", retries or 40))
    max_attempts = int(os.environ.get("GC_SEED_ATTEMPTS", max_attempts or 20000))
    reselect = bool(int(os.environ.get("GC_SEED_RESELECT",
                                       1 if reselect else 0)))

    # Default to IN-PROCESS seeding.
    #
    # The subprocess path below gives each attempt a hard wall-clock timeout,
    # which is the theoretically nicer design. In practice it did not work on
    # Windows: every attempt ran to the timeout and was killed, and six
    # instances burned 2 hours each producing nothing, even though the very
    # same settings seed FL in ~5 s in-process. Rather than debug Windows
    # multiprocessing, seed inline by default and keep the subprocess path
    # available with GC_SEED_SUBPROC=1.
    #
    # Inline seeding has no wall-clock timeout, so the attempt budget is kept
    # modest instead: a failure costs bounded work and we simply retry from a
    # new random seed, which is what actually helps given the heavy-tailed
    # runtime distribution.
    if os.environ.get("GC_SEED_SUBPROC", "0") != "1":
        return _seed_inline(json_path, k, pop_col, epsilon,
                            max_attempts, reselect, retries, log)

    ctx = mp.get_context("spawn")
    t_start = time.time()
    spawn_broken = False

    for attempt in range(retries):
        if spawn_broken:
            break
        seed = 12345 + 7919 * attempt
        q = ctx.Queue()
        p = ctx.Process(target=_worker,
                        args=(json_path, k, pop_col, epsilon,
                              max_attempts, reselect, seed, q))
        t0 = time.time()
        p.start()
        p.join(timeout)

        if p.is_alive():
            p.terminate()
            p.join()
            log("    seed attempt %d/%d: timed out after %ds, retrying"
                % (attempt + 1, retries, timeout))
            continue

        status = payload = None
        try:
            status, payload = q.get_nowait()
        except Exception:                                      # noqa: BLE001
            pass

        if status == "ok":
            log("    seeded on attempt %d/%d in %.1fs (total %.1fs)"
                % (attempt + 1, retries, time.time() - t0,
                   time.time() - t_start))
            return payload

        # On Windows, multiprocessing uses "spawn", which re-imports the calling
        # module. If the entry script is not guarded by
        #     if __name__ == "__main__":
        # the child dies instantly with a bootstrapping/freeze_support error and
        # every attempt would "fail" in a fraction of a second. Detect that and
        # fall back to seeding in-process, so a harness quirk can never look like
        # an infeasible instance.
        if payload is None and (time.time() - t0) < 5 and attempt >= 1:
            spawn_broken = True
            log("    subprocess seeding unavailable (spawn returned nothing); "
                "falling back to in-process seeding")
            break

        log("    seed attempt %d/%d failed after %.1fs (%s), retrying"
            % (attempt + 1, retries, time.time() - t0,
               payload or "no result returned"))

    if spawn_broken:
        return _seed_inline(json_path, k, pop_col, epsilon,
                            max_attempts, reselect, retries, log)

    log("    SEEDING FAILED after %d attempts (%.0f min total)"
        % (retries, (time.time() - t_start) / 60.0))
    return None


def _seed_inline(json_path, k, pop_col, epsilon, max_attempts, reselect,
                 retries, log):
    """Fallback path: seed in this process. No wall-clock timeout is possible
    here, so keep max_attempts modest and rely on retries from fresh seeds."""
    from functools import partial
    from gerrychain import Graph
    from gerrychain.tree import recursive_tree_part, bipartition_tree

    warnings.filterwarnings("ignore")
    G = Graph.from_json(json_path)
    target = sum(G.nodes[i][pop_col] for i in G.nodes()) / k
    method = partial(bipartition_tree, max_attempts=max_attempts,
                     allow_pair_reselection=reselect)
    for attempt in range(retries):
        t0 = time.time()
        try:
            random.seed(12345 + 7919 * attempt)
            a = recursive_tree_part(G, range(k), target, pop_col,
                                    epsilon, 1, method=method)
            log("    seeded in-process on attempt %d/%d in %.1fs"
                % (attempt + 1, retries, time.time() - t0))
            return {int(v): int(d) for v, d in a.items()}
        except Exception as exc:                               # noqa: BLE001
            log("    in-process seed attempt %d/%d failed after %.1fs (%s)"
                % (attempt + 1, retries, time.time() - t0, type(exc).__name__))
    log("    SEEDING FAILED (in-process) after %d attempts" % retries)
    return None


if __name__ == "__main__":
    # Smoke test:  python gc_seed.py CA 52
    import sys
    st, k = sys.argv[1], int(sys.argv[2])
    a = seed_partition("./raw_data/tract/json/%s_tracts.json" % st, k)
    print("RESULT:", "ok, %d districts" % len(set(a.values())) if a else "failed")
