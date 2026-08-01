# Smarter fixing ideas (for the paper and the code)

The current three-phase procedure decides each candidate vertex `v` by solving a
small feasibility MIP on the radius-`s` ego graph. Below are ideas that fix more
vertices, faster — one is already implemented and validated; the rest are
drop-in extensions worth a short subsection in the paper.

## 1. Combinatorial pre-screens — no MIP  *(IMPLEMENTED + validated)*

**Key fact.** If a district `D` has `diam(D) ≤ s` and `v ∈ D`, then every vertex
of `D` is within distance `s` of `v` inside `G[D]`, hence inside `G`. So
`D ⊆ B_s(v)`, the radius-`s` ball of `v`. This gives two *sound* necessary
conditions that need only a BFS and a knapsack — no solver.

> **Proposition (population ball).** If `pop(B_s(v)) < L`, then `v` belongs to no
> feasible district (diameter `≤ s`, population `≥ L`); fix `x_{v,·}=0`.

> **Proposition (minority surplus).** Let `surplus_i = p^{MVAP}_i − f·p^{VAP}_i`.
> If `max { Σ_{i∈S} surplus_i : v∈S⊆B_s(v), pop(S) ≥ L } < 0`, then no subset of
> `B_s(v)` — in particular no feasible district — is minority-majority; fix `v`.
> The maximum is taken over a relaxation (drop connectivity, the diameter, and
> the upper bound `U`), so its **LP value** is a sound over-estimate and is
> computed greedily in `O(|B_s(v)| log|B_s(v)|)`.

Implemented as `can_fix_minority(...)` in `fast_graph.py` and called as a
pre-screen in `fixing_sub_problem` (it fixes `v` and skips the MIP when the test
fires). Validated against brute force in `test_prefilter.py`: **0 false fixings**
over 1,300 vertices; it fires on ~6% of *random* instances and substantially
more on real rural / low-minority geographies, which is exactly where the MIP
solves were cheapest-to-avoid but most numerous. Because the pre-screen is a
relaxation of every pass, the set of fixings is unchanged — only faster.

**Paper angle:** these are clean, citable lemmas that make the fixing procedure
partly *certificate-based* (a vertex is fixed with a one-line combinatorial
proof, not a MIP run). Report the % of vertices resolved by the pre-screen vs the
MIP per state.

## 2. Monotonicity of fixings in `s`  *(easy proposition + reuse)*

> **Proposition.** A vertex fixed at diameter bound `s` is also fixable at every
> `s' ≤ s` (a smaller diameter is strictly more restrictive).

So when sweeping `s` (Section 5.2 / 6), fixings computed at one `s` transfer
downward for free — never recompute them. Likewise, fixings accumulate
monotonically across the fix-point iterations (removing fixed vertices only
shrinks each ego graph), which the code already exploits.

## 3. Reduced-cost fixing  *(TESTED — not worth it here)*

Two flavors, both ruled out:

- **In-model** reduced-cost fixing on the main MIP is something **Gurobi already
  does internally** (presolve + node bound strengthening), so re-implementing it
  adds nothing.
- **Cross-model**: take reduced costs from the relaxed upper-bound model R and
  fix variables in the full model P. This is *not* something Gurobi does (two
  different models), so it was worth a test. `test_crossmodel_fixing.py` runs it
  end to end and the result is decisive: on MS county (k=4, s=6) it fixed
  **0 of 328** assignment variables. Reasons: (i) R's LP bound is loose
  (`z_R_lp = 3.178` vs optimum `z* = 1`, gap ≈ 2.18), and (ii) the assignment
  variables `d[v,j]` have **zero objective coefficient**, so their reduced costs
  are ≈ 0 and never exceed the gap. The transfer is valid (the re-solve returns
  the same optimum) but empty.

  **Conclusion:** reduced-cost fixing is not a useful lever for this problem; the
  signal lives in the *combinatorial structure* (the radius-`s` ball), captured
  by the pre-screen in #1 — not in objective reduced costs. Worth one sentence in
  the paper to pre-empt the "why not reduced-cost fixing?" question.

## 4. Majority-side pre-screen  *(symmetric to #1)*

For the `option='majority'` pass, fix `v` when *every* pop-feasible subset of
`B_s(v)` containing `v` is minority-majority, i.e. when
`min { Σ surplus_i : v∈S, pop(S) ≥ L } > 0`. Same greedy LP with `min` instead of
`max`. (The current pipeline only runs the minority pass, so this is optional.)

## 5. One persistent model + warm starts  *(engineering)*

The per-vertex `nx.power`/`nx.complement` cost is already removed (`far_pairs`).
The remaining win is to stop rebuilding the Gurobi model per vertex:

1. one shared `gurobipy.Env()` passed to each subproblem; and
2. a single model with `t_u` for all `u`; per candidate set `t_v.LB=1` and
   `t_u.UB=0` outside `B_s(v)` (and for already-fixed `u`), re-solve, reset.
   The population/minority constraints are global sums that work unchanged; only
   the length-`s` conflict pairs are candidate-specific (add lazily).

Left unimplemented because it changes model construction (untestable without
Gurobi here); the pre-screen in #1 already removes most of the per-vertex solves,
so do this only if profiling still shows model build/solve overhead dominating.

## 6. Tighter knapsack (optional, closes a few more)

The pre-screen relaxation in #1 drops `U` and connectivity. Adding the upper
bound `U` (bounded knapsack) or a cheap connectivity surrogate (only count
vertices reachable from `v` within the budget) tightens it and fixes a few more
vertices, at higher cost. Worth trying only if the LP version leaves many
borderline vertices to the MIP.

---

### Suggested paper text (one subsection)
State Propositions in #1 and #2, report a table per state with columns:
*vertices, fixed by pre-screen (no MIP), fixed by MIP, total fixed (%), time*.
This directly answers the reviewers' interest in the fixing procedure's strength
while showing it is now largely certificate-based and much cheaper.
