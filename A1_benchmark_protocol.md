# A1 — Benchmark: our diameter model (full) vs. Belotti–Buchanan–Ezazipour published PP results

Addresses Reviewer 2, major comments 1, 2, and 5. We **re-run only our own
method** (the full optimized diameter model) on county-level and a handful of
small tract instances, score every reviewer-requested metric, and **compare
against the published results** of Belotti, Buchanan & Ezazipour (*Operations
Research* 2025), code/paper at
**https://github.com/AustinLBuchanan/Polsby_Popper_optimization** — we do **not**
re-run their code (they report their numbers; ours are what we add).

## Ground rules

- **Our side: the FULL optimized model** (all accelerations) — reduced
  formulation, three-phase parcel fixing, independent-set symmetry-breaking,
  up-front clique distance cuts, and the z-ordering. Objective: maximize the
  number of majority-minority districts subject to diameter ≤ s. We start at
  `s = ℓ_s` and increase until feasible (`s*`), with `TimeLimit = 3600 s`.
- **Their side: published numbers.** Take Belotti et al.'s reported compactness
  results (Polsby–Popper, and where available Reock / convex-hull) from their
  paper/repo and place them in the "Belotti et al. (published)" rows. Their method
  optimizes PP (it does not maximize #MM), so treat #MM there as reported/not
  applicable — note the difference in objective when discussing.
- **Same instances** (see the list below) — county level plus a handful of small
  tract instances. **Not** a full sweep.
- **Metrics** are computed on our returned plans the same way for every instance
  (see below): #MM, time, gap, Polsby–Popper, Reock, convex-hull, diameter.

## What to record (the metrics the reviewer asked for)

For each instance and each method, on the **final solution returned by that
method**, compute **all** of these with the **same** helpers (so the comparison
is apples-to-apples): (i) number of majority-minority districts, (ii) solve time
and final MIP gap, (iii) average **Polsby–Popper**, (iv) average **Reock**,
(v) average **convex-hull** score, and (vi) the **graph diameter** (max
within-district graph distance). Reviewer major comment 1 asks for #MM / time /
PP / diameter; major comment 2 asks for PP / Reock / convex-hull — all are
covered here and live in the workbook's columns.

The three geometry scores (PP, Reock, convex hull) are computed uniformly by
**dissolving the unit polygons into district polygons** (geopandas) and taking,
per district: `PP = 4π·area/perimeter²`, `convex-hull = area/area(convex hull)`,
`Reock = area/area(min bounding circle)` — see `benchmark_run_helpers.score_plan`
/ `geometry_scores`. (This geometric PP needs no `boundary_perim`.)

| State | level | n | k | method | #MM | time (s) | gap | avg PP | avg Reock | avg CH | max diam |
|-------|-------|---|---|--------|-----|----------|-----|--------|-----------|--------|----------|
| MS | county | 82 | 4 | diameter (ours, bare) | | | | | | | |
| MS | county | 82 | 4 | PP-MISOCP (theirs, bare) | | | | | | | |
| MS | county | 82 | 4 | PP-MISOCP (theirs, referee variant) | | | | | | | |
| ... | | | | | | | | | | | |

**Score every method's final plan identically.** For our runs this is automatic
(`run_A1_ours.py`). For their plans, export the assignment and re-score with the
same code:

```python
import benchmark_run_helpers as H
labels, districts, k = H.load_plan_csv("maps/MS_county_theirs_sol1.csv")
import geopandas
gdf = geopandas.read_file("raw_data/county/shape/MS_counties.shp")
print(H.score_plan(G, districts, mvap, vap, gdf=gdf, labels=labels))
# -> {'mm':…, 'diameter':…, 'avg_pp':…, 'avg_reock':…, 'avg_chull':…}
```

## Our side — FULL diameter model, "start at ℓ_s and increase until feasible"

Run the full optimized model with all accelerations (automated in
`run_A1_ours.py`):

```python
# inside the loop, after building `instance` for (state, level, group)
instance.model_time_limit = 3600
s = instance.lower_bound_s            # start at the lower bound on s
while True:
    instance.s = s
    instance.labeling_model(reduced_model=True,
                            parcel_fixing=True,
                            symmetry_breaking_constraint=True,
                            bare=False)                 # FULL: all accelerations
    if instance.num_minority_districts != 'Infeasible':
        break                                           # feasible at this s
    s += 1
```

`run_A1_ours.py` does exactly this for every instance, captures all feasible
solutions + their maps, scores all metrics, and writes the workbook.
(`parcel_fixing`/`symmetry_breaking_constraint` read precomputed fixing /
independent-set files when present and compute them otherwise — the latter can be
slow for tract instances that have not been preprocessed.)

## Their side — use Belotti et al.'s published results (no re-run)

We do **not** re-run their code. Take their reported numbers from the paper/repo
(`county_results/`) and place them in the "Belotti et al. (published)" rows of
the workbook. Two things to keep accurate when filling and discussing them:

- **Their objective differs.** Their method *optimizes* Polsby–Popper
  (their primary model minimizes the average inverse PP for a fixed `k`); it does
  **not** maximize the number of majority-minority districts. So #MM on their side
  is "as reported by their voting-rights study," not the result of an MM-maximizing
  objective — note this when comparing #MM.
- **Their reported plans use their full method**, including the Carve–Complete–
  Cleanup heuristic, warm starts, fixings, and orbitope symmetry — i.e. *their*
  best-tuned results, which is the right thing to compare our full model against.

If the published tables report Polsby–Popper but not Reock / convex-hull for a
given instance, leave those cells blank (or compute them only if they release the
plan). Cite the source table/figure next to the numbers.

> A direct head-to-head MISOCP run (and the referee's "their PP constraint + our
> MM objective" variant) is *not* part of this protocol since we are not re-running
> their code. We rely on their published compactness results instead, and on our
> own plans' PP / Reock / convex-hull scores vs. the enacted maps (major comment 2).

## Scoring any plan — one helper, all metrics (`score_plan`)

Score every method's final plan with the **same** function so the table is
symmetric: `benchmark_run_helpers.score_plan(G, districts, mvap, vap, gdf=gdf,
labels=labels)` returns `{mm, diameter, avg_pp, avg_reock, avg_chull}`.

- **# MM districts** and **diameter** need only networkx + `fast_graph` (no
  geopandas): diameter is `max_{a∈D} ecc(a)` via capped BFS in each district
  (`plan_diameter`), validated to match `nx.diameter` exactly.
- **Polsby–Popper, Reock, convex hull** are computed by `geometry_scores`, which
  dissolves the unit polygons into district polygons (geopandas) and takes
  `PP = 4π·area/perimeter²`, `convex-hull = area/area(convex hull)`,
  `Reock = area/area(min bounding circle)`. This geometric PP needs **no**
  `boundary_perim` (it uses the dissolved polygon's own area/perimeter), so it is
  simpler and more standard than the edge-based PP. Returns `None` gracefully if
  geopandas/shapely are unavailable.

## Instances (county + a handful of small tracts — not all)

- **County level** (small, solve fast): MS (82), AL (67), LA (64), SC, GA, MO —
  pick the ones already in the paper's county set.
- **A handful of small tract instances:** choose the few states with the *fewest*
  tracts in the paper's set (so the bare MISOCP is tractable), e.g. a couple of
  the smaller majority-Black tract states and one or two Latino tract states
  (NM, CO). Keep it to ~3–5 tract instances; this is a spot-check, not a sweep.

## Time limit + capturing ALL feasible solutions and their maps

- **One-hour time limit per MIP run:** set `TimeLimit = 3600` for every run
  (ours and theirs). If a run times out, record the incumbent objective and the
  **final MIP gap**.
- **Report every feasible solution found, not just the best.** Enable Gurobi's
  solution pool before solving:

  ```python
  m.Params.TimeLimit = 3600
  m.Params.PoolSearchMode = 2     # keep the n best solutions found
  m.Params.PoolSolutions  = 100   # retain up to 100 feasible solutions
  ```

  After `optimize`, iterate the pool and read each solution:

  ```python
  import fast_graph as fg
  sols = []
  for i in range(m.SolCount):
      m.Params.SolutionNumber = i
      plan = extract_districts(m, m._X, m._Y, n, k)   # use Xn for pool values
      # IMPORTANT (lazy constraints): pool entries are not guaranteed to satisfy
      # the lazy separator cuts, so VALIDATE each plan before reporting it.
      ok, why = fg.check_feasible_plan(csr, labels_of(plan), k, s, pop, L, U)
      if not ok:
          continue
      sols.append(plan)
  ```

  `fast_graph.check_feasible_plan` (tested) confirms each plan is a genuine
  feasible districting (partition, connected, diameter ≤ s, population in [L,U]).
  For the PP-MISOCP side, validate contiguity the same way (their plans need not
  satisfy a diameter bound — just report the diameter you measure).
- **Save each feasible plan as a map.** Write the node→district assignment to a
  CSV (and optionally render it from the shapefile), one file per solution, e.g.
  `maps/MS_county_ours_sol1.csv`. Record the filename in the workbook so every
  reported solution is traceable to its map.

## Output — one Excel workbook (`benchmark_to_excel.py`)

Use the provided `benchmark_to_excel.py` to write all metrics into a single,
formatted workbook with two sheets:

- **"A1 benchmark"** — one row per instance, best result of each method side by
  side: `# MM districts, time, gap, avg PP, max diameter`, plus `Δ # MM` and
  `time ratio` comparison columns.
- **"All feasible solutions"** — one row per feasible solution found within the
  hour (every validated pool entry), with its objective, time found, gap, PP,
  diameter, and the path to its saved map.

```python
from benchmark_to_excel import build_workbook
build_workbook(summary_rows, solution_rows, "A1_benchmark_results.xlsx")
```

Run `python benchmark_to_excel.py` with no arguments to emit a ready-to-fill
template (the instance list is pre-populated).

## Reading the result

The honest expected outcome (worth stating in the paper either way): the bare
PP-MISOCP is a continuous/SOC model Gurobi handles natively but it is typically
**harder** to solve to optimality than the bare diameter model at tract scale,
while the diameter plans should still score reasonably on PP. Whatever the
numbers show, this isolates the value of the compactness-measure choice exactly
as the referee requested.
