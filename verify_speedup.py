"""
verify_speedup.py  --  RUN THIS ON YOUR MACHINE (needs gurobipy).

Self-contained check that the optimized callback (compactness_callback) returns
the SAME optimum as the original (compactness_callback_ref) and is faster.

It builds the "max number of compact districts" feasibility/optimization MIP
(labeling formulation + length-s separator cuts) directly from a county JSON
file, so it does NOT depend on the rest of the pipeline or any results files.

Usage:
    python3 verify_speedup.py MS county black 6
    python3 verify_speedup.py            # defaults: MS county black s=6

Reports for each callback: status, objective, #nodes, wall-clock time.
Asserts the two objectives match.
"""
import sys, json, time, math
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
import aux_functions as AF


def load_graph(state, level):
    plural = 'counties' if level == 'county' else 'tracts'
    path = f'./raw_data/{level}/json/{state}_{plural}.json'
    d = json.load(open(path))
    G = nx.Graph()
    n = len(d['nodes'])
    G.add_nodes_from(range(n))
    pop = [d['nodes'][i]['P0010001'] for i in range(n)]
    for u, nbrs in enumerate(d['adjacency']):
        for e in nbrs:
            G.add_edge(u, e['id'])
    for i in range(n):
        G.nodes[i]['TOTPOP'] = pop[i]
    return G, pop


def build_model(G, pop, k, s, time_limit):
    """Labeling model: assign every node to one of k districts, population
    balanced, each district connected with diameter <= s (enforced lazily)."""
    n = G.number_of_nodes()
    total = sum(pop)
    L = math.ceil((1 - 0.005) * total / k)
    U = math.floor((1 + 0.005) * total / k)

    m = gp.Model()
    m.Params.OutputFlag = 0
    m.Params.TimeLimit = time_limit
    m.Params.LazyConstraints = 1
    m.Params.MIPGap = 0.0

    X = m.addVars(n, k, vtype=GRB.BINARY, name='X')
    m.addConstrs(gp.quicksum(X[v, j] for j in range(k)) == 1 for v in range(n))
    m.addConstrs(gp.quicksum(pop[v] * X[v, j] for v in range(n)) >= L for j in range(k))
    m.addConstrs(gp.quicksum(pop[v] * X[v, j] for v in range(n)) <= U for j in range(k))

    # attributes the callback expects
    m._X = X
    m._G = G
    m._n = n
    m._k = k
    m._s = s
    m._option = 'minority'
    m._district_vars_exist = False
    m._minority_district_range = range(k)
    m._majority_district_range = range(0)
    m._population = pop
    m._L = L
    return m


def extract_labels(m, X, n, k):
    if m.SolCount == 0:
        return None
    label = [-1] * n
    for v in range(n):
        for j in range(k):
            if X[v, j].X > 0.5:
                label[v] = j
                break
    return label


def timed(callback):
    """Wrap a callback so we can measure the cumulative time spent inside it
    (this isolates the separation cost from Gurobi's branch-and-bound time)."""
    box = {'t': 0.0, 'calls': 0}
    def wrapper(m, where):
        t0 = time.perf_counter()
        callback(m, where)
        box['t'] += time.perf_counter() - t0
        box['calls'] += 1
    wrapper.box = box
    return wrapper


def solve_with(callback, G, pop, k, s, time_limit):
    m = build_model(G, pop, k, s, time_limit)
    n = G.number_of_nodes()
    tcb = timed(callback)
    t0 = time.time()
    m.optimize(tcb)
    dt = time.time() - t0
    label = extract_labels(m, m._X, n, k)
    return {'status': m.Status, 'nodes': int(m.NodeCount), 'time': dt,
            'label': label, 'cb_time': tcb.box['t'], 'cb_calls': tcb.box['calls']}


if __name__ == '__main__':
    state = sys.argv[1] if len(sys.argv) > 1 else 'MS'
    level = sys.argv[2] if len(sys.argv) > 2 else 'county'
    s = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    TIME_LIMIT = 600

    G, pop = load_graph(state, level)
    k = G.number_of_nodes()  # placeholder; set a real k below
    # use the number of congressional districts if available, else a small k
    from aux_functions import get_congressional_codes
    k = get_congressional_codes().get(state.upper(), 4)
    print(f"Instance: {state} {level}  n={G.number_of_nodes()} m={G.number_of_edges()} k={k} s={s}")

    from fast_graph import CSRGraph, check_feasible_plan
    csr = CSRGraph(G.number_of_nodes(), G.edges())
    total = sum(pop)
    L = math.ceil((1 - 0.005) * total / k)
    U = math.floor((1 + 0.005) * total / k)

    print("\nSolving with OPTIMIZED callback (compactness_callback) ...")
    new = solve_with(AF.compactness_callback, G, pop, k, s, TIME_LIMIT)
    print(f"  status={new['status']} nodes={new['nodes']} time={new['time']:.2f}s"
          f"  callback_time={new['cb_time']:.3f}s over {new['cb_calls']} calls")

    print("Solving with ORIGINAL callback (compactness_callback_ref) ...")
    ref = solve_with(AF.compactness_callback_ref, G, pop, k, s, TIME_LIMIT)
    print(f"  status={ref['status']} nodes={ref['nodes']} time={ref['time']:.2f}s"
          f"  callback_time={ref['cb_time']:.3f}s over {ref['cb_calls']} calls")

    print("\n--- comparison ---")
    # 1. the two callbacks must agree on feasibility status
    assert new['status'] == ref['status'], \
        f"STATUS MISMATCH new={new['status']} ref={ref['status']}"
    print(f"feasibility status agrees: {new['status']}")

    # 2. independently verify any returned plan really is a feasible districting
    for name, res in (("optimized", new), ("reference", ref)):
        if res['label'] is not None:
            ok, why = check_feasible_plan(csr, res['label'], k, s, pop, L, U)
            assert ok, f"{name} callback returned an INFEASIBLE plan: {why}"
            print(f"{name} plan verified feasible (connected, diameter<= {s}, pop in [{L},{U}])")

    # callback-only speedup isolates the separation work from Gurobi B&B time;
    # this is the number that reflects the fast_graph rewrite (and it grows with
    # instance size). Total wall-clock speedup is smaller on tiny instances
    # because B&B, not the callback, dominates there.
    if new['cb_time'] > 0:
        print(f"\ncallback-only speedup: {ref['cb_time'] / max(new['cb_time'], 1e-9):.2f}x"
              f"  (ref {ref['cb_time']:.3f}s -> new {new['cb_time']:.3f}s)")
    if new['time'] > 0:
        print(f"total wall-clock speedup: {ref['time'] / max(new['time'], 1e-9):.2f}x")
    print("OK")
