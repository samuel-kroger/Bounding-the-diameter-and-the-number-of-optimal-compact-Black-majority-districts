"""
test_fast_graph_all.py  --  validate fast_graph.py against networkx.

Run from the code/ directory (needs only networkx, NOT gurobipy):

    python3 test_fast_graph_all.py

Checks, on random graphs, grids, and the real county adjacency graphs:
  1. CSR BFS distances == nx shortest path lengths
  2. components_within == nx connected components
  3. fischetti_separator == original find_fischetti_separator
  4. separate_compactness: every cut is a valid+minimal length-s separator,
     and exactly the violating districts are covered (so the optimized
     callback accepts the same integer solutions -> same optimum)
  5. separate_contiguity: valid full separators, all disconnected districts flagged
  6. power_edges == nx.power(G,s) edges; far_pairs == complement(power) pairs
"""
import os, json, random, glob
import networkx as nx
import fast_graph as fg


# ---------- reference helpers ---------- #
def ref_fischetti(DG, component, b, n):
    nb = [False] * n
    for i in nx.node_boundary(DG, component, None):
        nb[i] = True
    vis = [False] * n; child = [b]; vis[b] = True
    while child:
        par = child; child = []
        for i in par:
            if not nb[i]:
                for j in DG.neighbors(i):
                    if not vis[j]:
                        child.append(j); vis[j] = True
    return [i for i in DG.nodes if nb[i] and vis[i]]


def dist_minus(G, a, b, removed):
    H = G.copy(); H.remove_nodes_from([v for v in removed if v not in (a, b)])
    try:
        return nx.shortest_path_length(H, a, b)
    except nx.NetworkXNoPath:
        return float('inf')


def district_violates(G, Vj, s):
    if len(Vj) <= 1:
        return False
    H = G.subgraph(Vj)
    if not nx.is_connected(H):
        return True
    return nx.diameter(H) > s


# ---------- individual checks ---------- #
def check_bfs_components_fischetti(G):
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    DG = nx.to_directed(G)
    for src in list(G.nodes())[:min(15, n)]:
        g = csr.bfs1(src)
        ref = nx.single_source_shortest_path_length(G, src)
        for w in G.nodes():
            fd = csr._d1[w] if csr._s1[w] == g else None
            assert fd == ref.get(w, None)
    for _ in range(20):
        sub = set(random.sample(list(G.nodes()), random.randint(1, n)))
        member = [v in sub for v in range(n)]
        comps = fg.components_within(csr, sorted(sub), member)
        fast = sorted([frozenset(c) for c in comps], key=lambda z: (len(z), min(z)))
        ref = sorted([frozenset(c) for c in nx.connected_components(G.subgraph(sub))],
                     key=lambda z: (len(z), min(z)))
        assert fast == ref
        outside = [v for v in G.nodes() if v not in sub]
        if outside:
            b = random.choice(outside)
            for comp in nx.strongly_connected_components(DG.subgraph(sub)):
                cl = list(comp); mc = [v in comp for v in range(n)]
                assert sorted(fg._fischetti_for_component(csr, cl, b)) == \
                       sorted(ref_fischetti(DG, cl, b, n))


def check_separation(G, k, s, selector, trials):
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    mino = range(k) if selector != 'both' else range(max(1, k // 2))
    majo = range(0) if selector != 'both' else range(max(1, k // 2), k)
    for _ in range(trials):
        assign = [random.randrange(k) for _ in range(n)]
        xval = {(v, j): (1.0 if assign[v] == j else 0.0) for v in range(n) for j in range(k)}
        cuts = fg.separate_compactness(csr, s, selector, xval, xval, mino, majo, k)
        districts = {j: [v for v in range(n) if assign[v] == j] for j in range(k)}
        covered = set()
        for (a, b, minC) in cuts:
            assert assign[a] == assign[b]
            covered.add(assign[a])
            Vset = set(districts[assign[a]])
            assert all(c not in Vset for c in minC)
            assert dist_minus(G, a, b, set(minC)) > s
            for c in minC:
                assert dist_minus(G, a, b, set(minC) - {c}) <= s
        scan = mino if selector == 'minority' else range(k)
        for j in scan:
            if district_violates(G, districts[j], s):
                assert j in covered
        if not any(district_violates(G, districts[j], s) for j in scan):
            assert len(cuts) == 0


def check_contiguity(G, k, trials):
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    for _ in range(trials):
        assign = [random.randrange(k) for _ in range(n)]
        xval = {(v, j): (1.0 if assign[v] == j else 0.0) for v in range(n) for j in range(k)}
        cuts = fg.separate_contiguity(csr, 'minority', xval, {}, range(k), range(0), k)
        districts = {j: [v for v in range(n) if assign[v] == j] for j in range(k)}
        covered = set()
        for (a, b, C) in cuts:
            covered.add(assign[a])
            assert dist_minus(G, a, b, set(C)) == float('inf')
        for j in range(k):
            Vj = districts[j]
            if len(Vj) > 1 and not nx.is_connected(G.subgraph(Vj)):
                assert j in covered


def check_power(G):
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    sp = dict(nx.all_pairs_shortest_path_length(G))
    for s in [1, 2, 3]:
        ref = set(frozenset((u, v)) for u, v in nx.power(G, s).edges())
        fast = set(frozenset(e) for e in fg.power_edges(csr, s))
        assert ref == fast
        sub = sorted(random.sample(range(n), max(2, n // 2)))
        reff = set()
        for i, u in enumerate(sub):
            for v in sub[i + 1:]:
                d = sp[u].get(v, None)
                if d is None or d > s:
                    reff.add(frozenset((u, v)))
        assert reff == set(frozenset(e) for e in fg.far_pairs(csr, s, sub))


def load_county(path):
    d = json.load(open(path))
    G = nx.Graph(); G.add_nodes_from(range(len(d["nodes"])))
    for u, nbrs in enumerate(d["adjacency"]):
        for e in nbrs:
            G.add_edge(u, e["id"])
    return G


if __name__ == "__main__":
    random.seed(2024)
    print("random + grid graphs ...")
    graphs = []
    for i in range(15):
        g = nx.gnp_random_graph(random.randint(8, 40), random.uniform(0.12, 0.4), seed=i)
        if nx.is_connected(g):
            graphs.append(nx.convert_node_labels_to_integers(g))
    for (r, c) in [(5, 6), (7, 8)]:
        graphs.append(nx.convert_node_labels_to_integers(nx.grid_2d_graph(r, c)))
    for G in graphs:
        check_bfs_components_fischetti(G)
        check_power(G)
        for sel in ['minority', 'both']:
            for s in [2, 3]:
                check_separation(G, 4, s, sel, 5)
        check_contiguity(G, 4, 5)
    print("  OK")

    print("real county instances ...")
    for path in sorted(glob.glob("./raw_data/county/json/*_counties.json")):
        st = os.path.basename(path)[:2]
        G = load_county(path)
        if not nx.is_connected(G):
            continue
        check_bfs_components_fischetti(G)
        check_power(G)
        for sel in ['minority', 'both']:
            for s in [3, 5]:
                check_separation(G, 6, s, sel, 4)
        check_contiguity(G, 6, 4)
    print("  OK")
    print("\nALL FAST_GRAPH TESTS PASSED")
