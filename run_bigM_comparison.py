"""
run_bigM_comparison.py -- Reviewer 2, Minor comment 6.

Quantifies the advantage of the separate-variable majority-minority formulation
over the natural big-M alternative (in the spirit of Lawless and Gunluk 2024).
For each instance and each formulation we report the root LP-relaxation bound,
the number of branch-and-bound nodes, and the solve time.

To isolate the *formulation* effect, we compare on the core majority-minority
assignment problem -- maximize the number of majority-minority districts subject
to a partition, per-district population balance, the minority-share threshold, and
the z-ordering -- WITHOUT the diameter constraint (the diameter is enforced
identically in both formulations via the same length-s separator callback, so it
is not what differs).

  Separate-variable (ours):  x_{vj} for v, minority districts j  (z_j = 1 if MM)
                             y_{vj} for v, all districts j        (w_j = 1 if non-MM)
                             clean threshold  sum mvap*x >= f sum vap*x.
  Big-M (alternative):       a single x_{vj} for all districts, z_j in {0,1},
                             threshold enforced only when z_j = 1 via big-M:
                             sum mvap*x_{.j} - f sum vap*x_{.j} >= -M (1 - z_j).

The separate-variable subpolytope is integral (Proposition 2), so its LP
relaxation is tighter -- the table makes the gap, node-count, and time difference
explicit.

Run from code/ (needs gurobipy):  python run_bigM_comparison.py
Output: bigM_comparison.csv
"""
import csv
import json
import math
import time

import gurobipy as gp
from gurobipy import GRB

# small instances that solve quickly so both formulations finish
INSTANCES = [
    ("MS", "county", "counties", "black", 4),
    ("AL", "county", "counties", "black", 7),
    ("LA", "county", "counties", "black", 6),
]
DEVIATION = 0.01          # +/-0.5%, matching the paper
F = 0.5                   # minority-share threshold
TIME_LIMIT = 600

BLACK_ANYPART = [
    "P0030004", "P0030011", "P0030016", "P0030017", "P0030018", "P0030019",
    "P0030027", "P0030028", "P0030029", "P0030030", "P0030037", "P0030038",
    "P0030039", "P0030040", "P0030041", "P0030042", "P0030048", "P0030049",
    "P0030050", "P0030051", "P0030052", "P0030053", "P0030058", "P0030059",
    "P0030060", "P0030061", "P0030064", "P0030065", "P0030066", "P0030067",
    "P0030069", "P0030071",
]


def load(state, level, plural, group):
    with open("raw_data/%s/json/%s_%s.json" % (level, state, plural)) as fh:
        d = json.load(fh)
    nodes = d["nodes"]
    n = len(nodes)
    pop = [0] * n
    vap = [0] * n
    mvap = [0] * n
    for nd in nodes:
        i = nd["id"]
        pop[i] = nd["P0010001"]
        vap[i] = nd["P0030001"]
        mvap[i] = sum(nd[c] for c in BLACK_ANYPART) if group == "black" else nd["P0040002"]
    return n, pop, vap, mvap


def bounds(pop, k):
    ideal = sum(pop) / k
    return math.ceil((1 - DEVIATION / 2) * ideal), math.floor((1 + DEVIATION / 2) * ideal)


def build_separate(n, pop, vap, mvap, L, U, k):
    """Separate-variable formulation (ours)."""
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars(n, k, vtype=GRB.BINARY, name="x")   # node in minority district j
    z = m.addVars(k, vtype=GRB.BINARY, name="z")      # district j is MM
    y = m.addVars(n, k, vtype=GRB.BINARY, name="y")   # node in non-MM district j
    w = m.addVars(k, vtype=GRB.BINARY, name="w")      # district j is non-MM
    m.setObjective(gp.quicksum(z[j] for j in range(k)), GRB.MAXIMIZE)
    m.addConstrs(gp.quicksum(x[v, j] for j in range(k)) + gp.quicksum(y[v, j] for j in range(k)) == 1 for v in range(n))
    m.addConstrs(z[j] + w[j] == 1 for j in range(k))
    m.addConstrs(x[v, j] <= z[j] for v in range(n) for j in range(k))
    m.addConstrs(y[v, j] <= w[j] for v in range(n) for j in range(k))
    m.addConstrs(L * z[j] <= gp.quicksum(pop[v] * x[v, j] for v in range(n)) for j in range(k))
    m.addConstrs(gp.quicksum(pop[v] * x[v, j] for v in range(n)) <= U * z[j] for j in range(k))
    m.addConstrs(L * w[j] <= gp.quicksum(pop[v] * y[v, j] for v in range(n)) for j in range(k))
    m.addConstrs(gp.quicksum(pop[v] * y[v, j] for v in range(n)) <= U * w[j] for j in range(k))
    # clean minority-share threshold (no big-M):
    m.addConstrs(gp.quicksum(mvap[v] * x[v, j] for v in range(n)) >= F * gp.quicksum(vap[v] * x[v, j] for v in range(n)) for j in range(k))
    m.addConstrs(z[j] >= z[j + 1] for j in range(k - 1))   # z-ordering
    m.update()
    return m


def build_bigM(n, pop, vap, mvap, L, U, k):
    """Natural big-M alternative: a single assignment variable per district."""
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars(n, k, vtype=GRB.BINARY, name="x")   # node in district j (any)
    z = m.addVars(k, vtype=GRB.BINARY, name="z")      # district j is MM
    M = sum(vap)                                       # valid big-M (>= f * VAP_j)
    m.setObjective(gp.quicksum(z[j] for j in range(k)), GRB.MAXIMIZE)
    m.addConstrs(gp.quicksum(x[v, j] for j in range(k)) == 1 for v in range(n))
    m.addConstrs(L <= gp.quicksum(pop[v] * x[v, j] for v in range(n)) for j in range(k))
    m.addConstrs(gp.quicksum(pop[v] * x[v, j] for v in range(n)) <= U for j in range(k))
    # threshold enforced only when z_j = 1, via big-M:
    m.addConstrs(
        gp.quicksum(mvap[v] * x[v, j] for v in range(n)) - F * gp.quicksum(vap[v] * x[v, j] for v in range(n))
        >= -M * (1 - z[j]) for j in range(k))
    m.addConstrs(z[j] >= z[j + 1] for j in range(k - 1))   # z-ordering
    m.update()
    return m


def root_lp_bound(m):
    r = m.relax()
    r.Params.OutputFlag = 0
    r.optimize()
    return r.ObjVal if r.Status == GRB.OPTIMAL else None


def solve_ip(m):
    m.Params.OutputFlag = 0
    m.Params.TimeLimit = TIME_LIMIT
    m.Params.MIPGap = 0.0
    t0 = time.time()
    m.optimize()
    dt = time.time() - t0
    obj = m.ObjVal if m.SolCount > 0 else None
    return obj, int(m.NodeCount), dt, (m.MIPGap if m.SolCount > 0 else None)


def main():
    rows = []
    for (state, level, plural, group, k) in INSTANCES:
        print("\n=== %s %s (k=%d) ===" % (state, level, k))
        n, pop, vap, mvap = load(state, level, plural, group)
        L, U = bounds(pop, k)
        for name, builder in (("separate-variable (ours)", build_separate),
                              ("big-M", build_bigM)):
            m = builder(n, pop, vap, mvap, L, U, k)
            lp = root_lp_bound(m)
            opt, nodes, dt, gap = solve_ip(m)
            print("  %-26s  rootLP=%.3f  IP=%s  nodes=%d  time=%.1fs"
                  % (name, lp if lp is not None else float("nan"), opt, nodes, dt))
            rows.append(dict(state=state, level=level, n=n, k=k, formulation=name,
                             root_LP=round(lp, 3) if lp is not None else None,
                             IP_opt=opt, nodes=nodes, time_s=round(dt, 1),
                             gap=gap))
    with open("bigM_comparison.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=["state", "level", "n", "k", "formulation",
                                             "root_LP", "IP_opt", "nodes", "time_s", "gap"])
        wtr.writeheader()
        for r in rows:
            wtr.writerow(r)
    print("\nWrote bigM_comparison.csv")
    print("Expected: the separate-variable rootLP is closer to the integer optimum"
          " and solves in fewer nodes / less time than big-M.")


if __name__ == "__main__":
    main()
