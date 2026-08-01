"""
test_prefilter.py -- soundness of the no-MIP fixing pre-screen (fast_graph.can_fix_minority).

Run from code/ (needs only networkx):  python3 test_prefilter.py

For tiny instances it brute-forces all connected, pop-feasible, diameter-<=s
subsets containing v and checks: whenever can_fix_minority(v) returns True, NO
feasible minority-majority district contains v (i.e. the pre-screen never makes
a wrong fixing). Reports how often the pre-screen fires.
"""
import random, itertools
import networkx as nx
import fast_graph as fg


def feasible_district_exists(G, v, s, allowed, pop, L, U, minority, vap, f):
    H = G.subgraph([u for u in G.nodes() if allowed[u]])
    if v not in H:
        return False
    ball = list(nx.single_source_shortest_path_length(H, v, cutoff=s))
    if len(ball) > 18:
        return None  # skip (brute force too big)
    others = [u for u in ball if u != v]
    for r in range(len(others) + 1):
        for combo in itertools.combinations(others, r):
            S = (v,) + combo
            p = sum(pop[u] for u in S)
            if p < L or p > U:
                continue
            HS = H.subgraph(S)
            if not nx.is_connected(HS) or nx.diameter(HS) > s:
                continue
            if sum(minority[u] for u in S) >= f * sum(vap[u] for u in S):
                return True
    return False


def run(G, trials, seed):
    G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    rnd = random.Random(seed)
    checked = fired = violations = 0
    for _ in range(trials):
        pop = [rnd.randint(1, 4) for _ in range(n)]
        vap = [max(1, pop[i] - rnd.randint(0, 1)) for i in range(n)]
        minority = [rnd.randint(0, vap[i]) for i in range(n)]
        f = 0.5
        tot = sum(pop); k = rnd.randint(2, 4)
        L = max(1, tot // k - 2); U = tot // k + 2
        allowed = [True] * n
        for v in range(n):
            s = rnd.choice([2, 3])
            ref = feasible_district_exists(G, v, s, allowed, pop, L, U, minority, vap, f)
            if ref is None:
                continue
            fix = fg.can_fix_minority(csr, v, s, allowed, pop, L, U, minority, vap, f)
            checked += 1
            if fix:
                fired += 1
            if fix and ref is True:
                violations += 1
    return checked, fired, violations


if __name__ == "__main__":
    random.seed(5)
    tc = tf = tv = 0
    for i in range(40):
        G = nx.gnp_random_graph(random.randint(7, 13), random.uniform(0.2, 0.45), seed=i)
        if nx.is_connected(G):
            c, f, v = run(G, 4, seed=i); tc += c; tf += f; tv += v
    for (r, cc) in [(3, 4), (4, 4)]:
        c, f, v = run(nx.grid_2d_graph(r, cc), 6, seed=99); tc += c; tf += f; tv += v
    print("checked=%d  pre-screen fixed=%d (%.1f%%)  soundness violations=%d"
          % (tc, tf, 100 * tf / max(tc, 1), tv))
    assert tv == 0, "UNSOUND: pre-screen fixed a vertex that has a feasible district!"
    print("PRE-SCREEN IS SOUND (no false fixings)")
