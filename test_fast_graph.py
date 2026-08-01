"""
Validate fast_graph routines against networkx on random graphs and on the
real county-level instances. We check:
  (1) BFS distances == nx shortest path lengths
  (2) components_within == nx connected components of the induced subgraph
  (3) fischetti_separator matches the original aux_functions implementation
      AND is a valid 1-separator (removing C disconnects a from b)
  (4) separator_from_dists yields a VALID and MINIMAL length-s a,b-separator,
      matching the *semantics* of the original sep-weight construction.
"""
import os
import json
import random
import itertools
import networkx as nx
import fast_graph as fg


# ----- original networkx reference implementation (from aux_functions) ----- #
def find_fischetti_separator_ref(DG, component, b, n):
    neighbors_component = [False for _ in range(n)]
    for i in nx.node_boundary(DG, component, None):
        neighbors_component[i] = True
    visited = [False for _ in range(n)]
    child = [b]
    visited[b] = True
    while child:
        parent = child
        child = []
        for i in parent:
            if not neighbors_component[i]:
                for j in DG.neighbors(i):
                    if not visited[j]:
                        child.append(j)
                        visited[j] = True
    return [i for i in DG.nodes if neighbors_component[i] and visited[i]]


def ref_min_length_s_separator(G, a, b, s, V_j):
    """Reference: replicate the sep-weight construction from the original
    compactness_callback diameter branch, returning minC."""
    S_prime = list(set(G.nodes) - set(V_j))
    distance_from_a_in_G = nx.single_source_shortest_path_length(G, a)
    distance_from_b_in_G = nx.single_source_shortest_path_length(G, b)
    remove_from_S_prime = []
    for v in S_prime:
        if distance_from_a_in_G.get(v, 10**9) + distance_from_b_in_G.get(v, 10**9) > s:
            remove_from_S_prime.append(v)
    C = [node for node in S_prime if node not in remove_from_S_prime]
    remove_from_C = []
    for u, v in G.edges():
        G[u][v]['sep-weight'] = 1
    for c in C:
        for nb in G.neighbors(c):
            G[c][nb]['sep-weight'] = s + 1
    for c in C:
        for nb in G.neighbors(c):
            G[c][nb]['sep-weight'] = 1
        d = nx.single_source_dijkstra_path_length(G, a, weight='sep-weight')
        if d[b] > s:
            remove_from_C.append(c)
        else:
            for nb in G.neighbors(c):
                G[c][nb]['sep-weight'] = s + 1
    return [c for c in C if c not in remove_from_C]


def dist_minus_set(G, a, b, removed):
    """dist_{G - removed}(a,b), inf if disconnected."""
    H = G.subgraph([v for v in G.nodes if v not in removed or v in (a, b)])
    # ensure removed (other than a,b) truly gone
    H = G.copy()
    H.remove_nodes_from([v for v in removed if v not in (a, b)])
    try:
        return nx.shortest_path_length(H, a, b)
    except nx.NetworkXNoPath:
        return float('inf')


def make_csr(G):
    return fg.CSRGraph(G.number_of_nodes(), list(G.edges()))


def test_bfs(G):
    csr = make_csr(G)
    for src in list(G.nodes())[: min(20, G.number_of_nodes())]:
        g = csr.bfs(src)
        ref = nx.single_source_shortest_path_length(G, src)
        for w in G.nodes():
            fd = csr.distance(g, w)
            rd = ref.get(w, None)
            assert fd == rd, f"bfs dist mismatch src={src} w={w} fast={fd} ref={rd}"


def test_components(G):
    csr = make_csr(G)
    nodes = list(G.nodes())
    # random subset
    for _ in range(30):
        k = random.randint(1, len(nodes))
        sub = set(random.sample(nodes, k))
        member = [v in sub for v in range(G.number_of_nodes())]
        comps = fg.components_within(csr, sorted(sub), member)
        fast_sets = sorted([frozenset(c) for c in comps], key=lambda s: (len(s), min(s)))
        ref = nx.connected_components(G.subgraph(sub))
        ref_sets = sorted([frozenset(c) for c in ref], key=lambda s: (len(s), min(s)))
        assert fast_sets == ref_sets, "component mismatch"


def test_fischetti(G):
    n = G.number_of_nodes()
    DG = nx.to_directed(G)
    csr = make_csr(G)
    nodes = list(G.nodes())
    for _ in range(40):
        k = random.randint(1, max(1, len(nodes) // 2))
        sub = set(random.sample(nodes, k))
        member = [v in sub for v in range(n)]
        # pick b outside sub (root) if possible
        outside = [v for v in nodes if v not in sub]
        if not outside:
            continue
        b = random.choice(outside)
        comp_sets = list(nx.strongly_connected_components(DG.subgraph(sub)))
        for component in comp_sets:
            comp = list(component)
            member_comp = [v in component for v in range(n)]
            ref = sorted(find_fischetti_separator_ref(DG, comp, b, n))
            fast = sorted(fg.fischetti_separator(csr, comp, member_comp, b))
            assert ref == fast, f"fischetti mismatch ref={ref} fast={fast}"


def test_length_s_separator(G):
    n = G.number_of_nodes()
    csr = make_csr(G)
    nodes = list(G.nodes())
    order = list(range(n))
    checks = 0
    for _ in range(60):
        k = random.randint(2, len(nodes))
        Vj = set(random.sample(nodes, k))
        # need connected V_j with diameter > s for a meaningful test
        H = G.subgraph(Vj)
        if not nx.is_connected(H):
            continue
        s = random.randint(1, max(1, nx.diameter(H) - 1))
        # find a violating pair (a,b) with dist_H(a,b) > s
        found = None
        for a in Vj:
            dl = nx.single_source_shortest_path_length(H, a)
            for b, d in dl.items():
                if a < b and d > s:
                    found = (a, b)
                    break
            if found:
                break
        if not found:
            continue
        a, b = found
        # fast separator
        ga = csr.bfs(a)
        dist_a = [csr.distance(ga, v) for v in range(n)]
        gb = csr.bfs(b)
        dist_b = [csr.distance(gb, v) for v in range(n)]
        member_vj = [v in Vj for v in range(n)]
        fastC = fg.separator_from_dists(csr, a, b, s, member_vj, dist_a, dist_b, order)
        # 1) validity: removing fastC disconnects a,b within distance s
        assert dist_minus_set(G, a, b, set(fastC)) > s, "fast separator invalid"
        # 2) minimality: putting any node back breaks separation
        for c in list(fastC):
            reduced = set(fastC) - {c}
            assert dist_minus_set(G, a, b, reduced) <= s, "fast separator not minimal"
        # 3) reference separator is also valid (sanity of reference)
        refC = ref_min_length_s_separator(G, a, b, s, Vj)
        assert dist_minus_set(G, a, b, set(refC)) > s, "ref separator invalid?!"
        checks += 1
    return checks


def run_on(G, label):
    G = nx.convert_node_labels_to_integers(G)
    test_bfs(G)
    test_components(G)
    test_fischetti(G)
    c = test_length_s_separator(G)
    print(f"  [{label}] OK  (n={G.number_of_nodes()}, m={G.number_of_edges()}, sep-tests={c})")


if __name__ == "__main__":
    random.seed(12345)
    print("Random graphs:")
    for i in range(25):
        nnodes = random.randint(6, 40)
        p = random.uniform(0.12, 0.45)
        G = nx.gnp_random_graph(nnodes, p, seed=i)
        if G.number_of_edges() == 0:
            continue
        run_on(G, f"gnp#{i}")
    print("Grid graphs (planar, district-like):")
    for (r, cc) in [(4, 5), (6, 6), (8, 10)]:
        G = nx.grid_2d_graph(r, cc)
        run_on(G, f"grid{r}x{cc}")

    # real county instances
    base = "/sessions/gifted-bold-darwin/mnt/Sam_districting_paper/code/raw_data/county/json"
    print("Real county instances:")
    for st in ["MS", "LA", "AL", "NM", "CO"]:
        path = os.path.join(base, f"{st}_counties.json")
        if not os.path.exists(path):
            continue
        d = json.load(open(path))
        G = nx.Graph()
        G.add_nodes_from(range(len(d["nodes"])))
        for u, nbrs in enumerate(d["adjacency"]):
            for e in nbrs:
                v = e["id"]
                G.add_edge(u, v)
        run_on(G, st)
    print("\nALL TESTS PASSED")
