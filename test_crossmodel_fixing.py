"""
test_crossmodel_fixing.py  --  RUN ON YOUR MACHINE (needs gurobipy).

Empirical test of CROSS-MODEL reduced-cost fixing:

  R = relaxed "max number of minority-majority districts" model WITHOUT
      contiguity/diameter (this is the Section-5.1-style upper bound).
  P = the full model = R + contiguity + diameter (length-s separator callback).

We solve the LP relaxation of R, read its reduced costs and LP objective z_R_lp
(a valid upper bound for P, because P is more constrained than R). Given a valid
lower bound z_LB for P (here we use the true optimum of P, the most generous
incumbent, to measure the maximum fixing power), a binary variable can be fixed
at its LP bound if   z_R_lp - |reduced_cost| < z_LB.

Gurobi does NOT do this automatically because the reduced costs come from a
DIFFERENT model (R) than the one being solved (P).

The script then re-solves P with those fixings and ASSERTS the optimum is
unchanged, and reports how many assignment variables were fixed and the effect
on solve time / B&B nodes.

Usage:
    python3 test_crossmodel_fixing.py MS county black 6
    python3 test_crossmodel_fixing.py AL county black 6
"""
import sys, json, math, time
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
import aux_functions as AF
from aux_functions import get_congressional_codes


# ----------------------------- data ----------------------------- #
def load(state, level, group):
    plural = 'counties' if level == 'county' else 'tracts'
    d = json.load(open(f'./raw_data/{level}/json/{state}_{plural}.json'))
    n = len(d['nodes'])
    G = nx.Graph(); G.add_nodes_from(range(n))
    for u, nbrs in enumerate(d['adjacency']):
        for e in nbrs:
            G.add_edge(u, e['id'])
    pop = [d['nodes'][i]['P0010001'] for i in range(n)]
    vap = [d['nodes'][i]['P0030001'] for i in range(n)]
    code = 'P0030004' if group == 'black' else 'P0040002'
    mvap = [d['nodes'][i][code] for i in range(n)]
    return G, pop, vap, mvap


# ----------------- build the max-minority model ----------------- #
def build(G, pop, vap, mvap, k, s, f, contiguity, relax=False):
    """Separate-variable (big-M-free) formulation:
       d[v,j] = v assigned to district j ; x = d if j minority-majority else 0.
       max sum_j z[j].  If contiguity=True, length-s separator cuts on d."""
    n = G.number_of_nodes()
    total = sum(pop)
    L = math.ceil((1 - 0.005) * total / k)
    U = math.floor((1 + 0.005) * total / k)
    vt = GRB.CONTINUOUS if relax else GRB.BINARY

    m = gp.Model(); m.Params.OutputFlag = 0
    d = m.addVars(n, k, vtype=vt, lb=0, ub=1, name='d')
    x = m.addVars(n, k, lb=0, ub=1, name='x')      # minority part of d
    z = m.addVars(k, vtype=vt, lb=0, ub=1, name='z')

    m.addConstrs(gp.quicksum(d[v, j] for j in range(k)) == 1 for v in range(n))
    m.addConstrs(gp.quicksum(pop[v] * d[v, j] for v in range(n)) >= L for j in range(k))
    m.addConstrs(gp.quicksum(pop[v] * d[v, j] for v in range(n)) <= U for j in range(k))
    # x = d when z=1, x = 0 when z=0  (x <= d, x <= z, x >= d - (1-z))
    m.addConstrs(x[v, j] <= d[v, j] for v in range(n) for j in range(k))
    m.addConstrs(x[v, j] <= z[j] for v in range(n) for j in range(k))
    m.addConstrs(x[v, j] >= d[v, j] - (1 - z[j]) for v in range(n) for j in range(k))
    # minority-majority on the x (minority) part of each district
    m.addConstrs(gp.quicksum((mvap[v] - f * vap[v]) * x[v, j] for v in range(n)) >= 0
                 for j in range(k))

    m.setObjective(gp.quicksum(z[j] for j in range(k)), GRB.MAXIMIZE)
    m._d = d; m._x = x; m._z = z; m._L = L; m._U = U

    if contiguity:
        m.Params.LazyConstraints = 1
        m._X = d                       # callback enforces contiguity+diameter on d
        m._G = G; m._n = n; m._k = k; m._s = s
        m._option = 'minority'; m._district_vars_exist = False
        m._minority_district_range = range(k); m._majority_district_range = range(0)
    return m, L, U


def solve_P(G, pop, vap, mvap, k, s, f, time_limit, fix0=None):
    m, L, U = build(G, pop, vap, mvap, k, s, f, contiguity=True, relax=False)
    m.Params.TimeLimit = time_limit
    m.Params.MIPGap = 0.0
    if fix0:
        for (v, j) in fix0:
            m._d[v, j].ub = 0.0
    t0 = time.time()
    m.optimize(AF.compactness_callback)
    return m, time.time() - t0


if __name__ == '__main__':
    state = sys.argv[1] if len(sys.argv) > 1 else 'MS'
    level = sys.argv[2] if len(sys.argv) > 2 else 'county'
    group = sys.argv[3] if len(sys.argv) > 3 else 'black'
    s = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    f = 0.5
    TL = 1200

    G, pop, vap, mvap = load(state, level, group)
    n = G.number_of_nodes()
    k = get_congressional_codes().get(state.upper(), 4)
    print(f"{state} {level} {group}: n={n} m={G.number_of_edges()} k={k} s={s}")

    # 1) LP relaxation of the RELAXED model R (no contiguity) -> reduced costs
    R, _, _ = build(G, pop, vap, mvap, k, s, f, contiguity=False, relax=True)
    R.optimize()
    z_R_lp = R.ObjVal
    rc = {(v, j): R._d[v, j].RC for v in range(n) for j in range(k)}
    dval = {(v, j): R._d[v, j].X for v in range(n) for j in range(k)}
    print(f"R LP relaxation upper bound z_R_lp = {z_R_lp:.3f}")

    # 2) true optimum of P (valid lower bound; most generous incumbent)
    P0, t0 = solve_P(G, pop, vap, mvap, k, s, f, TL)
    if P0.SolCount == 0:
        print("P found no solution in the time limit; pick an easier instance.")
        sys.exit(0)
    z_star = round(P0.ObjVal)
    print(f"P optimum z* = {z_star}  (status {P0.Status}, {t0:.1f}s, {int(P0.NodeCount)} nodes)")

    # 3) cross-model reduced-cost fixing of d[v,j] at their LP bound 0
    fix0 = []
    for (v, j), val in dval.items():
        if val < 1e-6:                       # d[v,j] nonbasic at LB 0
            if z_R_lp - abs(rc[(v, j)]) < z_star - 1e-6:
                fix0.append((v, j))
    tot = n * k
    print(f"cross-model RC fixing: {len(fix0)} of {tot} assignment vars fixed to 0"
          f"  ({100*len(fix0)/tot:.1f}%)")

    if not fix0:
        print("=> no variables fixable on this instance (expected when R's LP is loose"
              " or assignment reduced costs are ~0).")
        sys.exit(0)

    # 4) re-solve P with the fixings; verify same optimum, compare effort
    P1, t1 = solve_P(G, pop, vap, mvap, k, s, f, TL, fix0=fix0)
    z1 = round(P1.ObjVal) if P1.SolCount else None
    print(f"P with fixings: z = {z1}  ({t1:.1f}s, {int(P1.NodeCount)} nodes)")
    assert z1 == z_star, "OPTIMUM CHANGED -- cross-model fixing was INVALID!"
    print("optimum unchanged (fixing valid).")
    print(f"time: {t0:.1f}s -> {t1:.1f}s   nodes: {int(P0.NodeCount)} -> {int(P1.NodeCount)}")
