# Code speedups — callback, power graphs, fixing, planarity

This summarizes the performance changes, how they were validated, and the
planarity analysis. **All changes preserve the model and the optimum** — the
optimized separation accepts exactly the same integer solutions as the original.

## What changed

### 1. The separation callback (the main hot spot) — `aux_functions.py`
`compactness_callback` and `continuous_ONLY_callback` were rewritten on top of a
new array-based module, **`fast_graph.py`**. The originals are kept verbatim as
`compactness_callback_ref` / `continuous_ONLY_callback_ref` (drop-in fallbacks
and the reference for testing).

The original callback, on **every** `MIPSOL`, rebuilt the directed graph
(`nx.to_directed`), and per district called `nx.is_strongly_connected`,
`nx.diameter` (all-pairs), and repeated `nx.single_source_dijkstra_path_length`
with the `'sep-weight'` edge trick to minimalize separators. That is where the
time went.

`fast_graph.py` replaces all of this with CSR adjacency + array BFS:
- the CSR is **built once** and cached on the model (`m._csr`), not rebuilt per call;
- one O(n·k) pass reads the solution; connectivity uses a BFS component scan;
- diameter violations are found by an **early-exit** BFS (no full all-pairs diameter);
- minimal length-s separators are built with a capped BFS + greedy reduction,
  using three independent scratch areas (generation-stamp trick) so nothing is
  reallocated between separations;
- the disconnected branch uses a fast Fischetti separator (same as the original),
  then length-s minimalization.

The "sep-weight Dijkstra > s" test is provably equivalent to "BFS distance in
G with the separator nodes removed > s", which is what the fast code does.

### 2. Power-graph construction — `aux_functions.compute_stability_number`, `comp_experiment.add_lazy_distance_cuts`, `...add_lazy_distance_cuts_fixing_problem`
`nx.power(G, s)` (and `nx.complement(...)`) materialize a dense graph. They are
replaced by two BFS helpers in `fast_graph.py`:
- `power_edges(csr, s)` — edges of G^s (pairs within distance s);
- `far_pairs(csr, s, nodes)` — pairs with distance > s (= complement of G^s).

Both are validated **equal** to the networkx versions.

### 3. Fixing subproblem — `comp_experiment.add_lazy_distance_cuts_fixing_problem`
Now uses `far_pairs` instead of `nx.power`+`nx.complement` per candidate vertex.

## Benchmarks (separation only, full cut construction on both sides)

Identical "flagged district" sets in every trial (so identical accept/reject ->
identical optimum). Random + realistic near-feasible solutions:

| instance | n | speedup (new vs original networkx) |
|---|---|---|
| MS / AL / LA / CO counties | 64–82 | **6–8×** |
| grid 30×30 | 900 | **18×** |
| grid 45×45 | 2025 | **21×** |
| grid 60×60 (tract scale) | 3600 | **24×** |

The gap **widens with size** because the original cost is dominated by
`nx.diameter` / Dijkstra, which scale poorly, while the new code is near-linear
per separation.

## How it was validated (no Gurobi needed for these)

Run from `code/`:

```
python3 test_fast_graph_all.py
```

This checks, against networkx, on random graphs, grids, and every real county
graph:
1. BFS distances == nx shortest paths;
2. connected components == nx;
3. Fischetti separator == the original implementation;
4. **every** cut from `separate_compactness` is a valid **and minimal** length-s
   a,b-separator, and **exactly** the violating districts are flagged (so no
   infeasible solution is ever accepted, and feasible solutions yield no cuts);
5. contiguity separators are valid and complete;
6. `power_edges` / `far_pairs` == `nx.power` / `nx.complement`.

The Gurobi wrapper logic (cut emission) was checked with a mock model that
records `cbLazy` calls.

### End-to-end check on your machine (needs gurobipy)

```
python3 verify_speedup.py MS county black 6
```

It builds the labeling + length-s-separator MIP directly from a county JSON,
solves it once with the **new** callback and once with the **reference**
callback, and **asserts the objectives match**, reporting both wall-clock times.
Please run this on a couple of county instances as the final gate before trusting
the new callback on the big tract runs.

## Planarity analysis (where it does / doesn't help)

County/tract adjacency graphs are **planar**, so m ≤ 3n−6 = O(n). Implications:

- **BFS / separation are O(n).** The CSR rewrite already turns every BFS and
  separator search into a linear-time array scan. This is the main planarity win
  and it is realized.
- **Full diameter is avoided.** Planar diameter has specialized algorithms, but
  the callback never needs the diameter value — it early-exits on the first
  violating pair, which is strictly cheaper.
- **Power graph G^s is *not* sparse even when G is planar.** A vertex's
  s-neighborhood can cover most of the graph, so G^s and its complement can have
  Θ(n²) edges. Planarity does **not** shrink them. What helps is building them by
  BFS (`power_edges` / `far_pairs`) instead of `nx.power`, avoiding the dense
  intermediate object and networkx overhead — done.
- **Larger, optional lever (not implemented):** Lipton–Tarjan planar separators
  give O(√n)-size balanced separators, which could drive a divide-and-conquer for
  the lower-bound-s independent-set computation or a planar-aware warm start.
  This is a research-grade change; it should be prototyped and benchmarked
  separately rather than dropped into the solve path. Flagging it as the main
  remaining planarity opportunity.

## Fixing code — further opportunity (needs your Gurobi machine to validate)

You suggested "build one model and run the fixing code multiple times." The
graph-side cost (the per-vertex `nx.power`/`nx.complement`) is now removed via
`far_pairs`. The remaining idea — a single **persistent** Gurobi model reused
across candidate vertices — is sound but changes model construction, which I
could not test here (no Gurobi in this environment). Recommended safe version:

1. Create one `gurobipy.Env()` and pass `env=` to each subproblem model (removes
   repeated environment setup).
2. Keep one model with `t[v]` for all v; per candidate, set `t[center].LB = 1`
   and `t[u].UB = 0` for u outside the radius-s ego graph (and for already-fixed
   u), then re-solve; reset bounds afterward. The population / minority
   constraints are global sums and work unchanged; the only per-candidate cuts
   are the length-s conflict pairs, which can be added as lazy and left in place.

I left the model structure unchanged so nothing is broken; this refactor is the
next step once you can run the Gurobi regression (`verify_speedup.py` pattern).

## Note on files

`fast_graph.py`, `aux_functions.py`, `comp_experiment.py`, `verify_speedup.py`,
and `test_fast_graph_all.py` in this folder are the up-to-date versions. (You
have your own backup as requested.)
