# Reproducing the paper's tables (final code)

This directory regenerates every table in the paper with the **final** code, one
table per batch file, so the whole suite reproduces on a single machine. The old
workflow (hand-editing `comp_experiment.py`'s `__main__` to set `group`/`dataset`
and uncomment a method) is replaced by a parameterized driver, `run_table.py`.

## 0. Environment (once)

From an **Anaconda Prompt**:

```
conda activate base
pip install gurobipy gerrychain geopandas tqdm
```

You also need a working **Gurobi license** (academic is fine). Everything is run
from this `code\` directory.

## 1. Run the tables, in this order

Each `.bat` activates `base`, `cd`s here, runs both the Black and Latino instance
sets, and writes a CSV under `results_black\` / `results_hispanic\`. Run them in
order — later tables reuse the bounds cached by the earlier ones.

| Order | Batch file | Produces | CSV |
|------|------------|----------|-----|
| 1 | `table1_instance_info.bat` | instance sizes, diameter, % minority VAP | `instance_info.csv` |
| 2 | `table2_lower_bound_s.bat` | lower bound on `s` (`ell_s`) | `lower_bound_s_table.csv` |
| 3 | `table3_upper_bound_mm.bat` | upper bound on #majority-minority districts + % fixed | `upper_bound_minority_table.csv` |
| 4 | `table4_fixing.bat` | three-phase fixing: lower-bound gains + % fixed per phase | `fixing_info.csv` |
| 5 | `table5_symmetry.bat` | symmetry-breaking speedup (with vs. without) | `symmetry_table.csv` |
| 6 | `table6_tradeoff.bat` | `(s, #MM)` tradeoff frontier (R1 C3 / R2 Major 3) | `point_check.csv` |
| 7 | `table7_benchmark.bat` | comparison vs. Belotti (PP-optimal, CCC) + GerryChain | `A1_benchmark_combined.xlsx` |
| 8 | `table8_bigM.bat` | separate-variable vs. big-M formulation (root LP, nodes, time) | `bigM_comparison.csv` |

Run 2 before 3/4/5: the instance constructor caches `ell_s` and the
independent-set files, which those tables reuse. Run 4 before 5 for the same
reason.

### The `run_table.py` driver (used by every batch)

```
python run_table.py --table <name> --group <black|hispanic|none> --dataset <key>
```
`--dataset` is a key in `data.json` (`paper` = the 17 majority-Black instances,
`hispanic_paper` = the 7 majority-Latino instances, `county` = the county-level
set). Run `python run_table.py -h` for the table list. The driver loops the
instances, builds each `problem_instance`, calls the matching method, and prints a
clean per-instance log; an error on one instance is caught and reported without
killing the run.

## 2. Table 7 (benchmark) — generate the competitors' plans first

`table7_benchmark.bat` re-runs **our** exact model and rebuilds the combined
workbook, but it needs the competitors' plans to exist. Generate them once:

```
cd ..\Polsby_Popper_optimization-main
python build_belotti_data.py          REM Belotti-format data from our graphs
run_belotti_our_instances.bat         REM Belotti PP-optimal MISOCP (county)
run_belotti_tracts.bat                REM Belotti PP-optimal MISOCP (tracts)
run_gingles.bat                       REM Belotti Carve-Complete-Cleanup (MM heuristic)
cd ..\code
python run_gerrychain.py              REM GerryChain upper bounds + warm starts
```
Then run `table7_benchmark.bat`. (The Belotti and GerryChain plan files are
deterministic given the data, so they need to be generated only once; only our
exact model must be re-run with the final code.)

## 3. Table 8 (big-`M` comparison) -- standalone

`table8_bigM.bat` runs `run_bigM_comparison.py`, which builds both the
separate-variable formulation and the natural big-`M` alternative (single `x`
plus a big-`M` minority-share constraint) on a few small instances and reports
the root LP-relaxation bound, B&B node count, and solve time for each. To isolate
the formulation effect it compares on the core majority-minority assignment
problem (no diameter constraint -- the diameter is enforced identically in both
via the same separator callback). It is standalone (not part of `run_table.py`).

## 4. Efficiency notes (what changed this round)

The repeated, per-call operations are the bottleneck, and they now use the
array-based routines in `fast_graph.py` instead of networkx:

* the length-`s` separator **callback** (contiguity + diameter) — `separate_*`;
* the **three-phase fixing** radius-`s` ball — `bfs1(..., cap=s)` (was
  `nx.ego_graph` per node) and a combinatorial pre-screen (`can_fix_minority`)
  that avoids a MIP for most parcels;
* the **independence-number / lower-`s`** certificate — `power_edges` (was
  `nx.power` + `nx.complement`);
* the up-front **distance cuts** — `far_pairs` (was `nx.power`), with the majority
  cuts measured on the full graph (a correctness fix this round);
* validated 4–24x speedups on these hot paths.

We deliberately **kept `nx.diameter`** for the one-time graph-diameter computation:
we implemented and tested an array version (`fast_graph.graph_diameter`, exact and
validated against `nx.diameter`), but networkx's implementation is actually faster
for a single diameter call, so it remains the right choice there.
