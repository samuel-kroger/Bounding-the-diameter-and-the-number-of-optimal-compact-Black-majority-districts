"""
Test separate_compactness() at the level of *cuts produced*:

  - SOUNDNESS: every produced (a,b,minC) has a,b in the same district, minC
    disjoint from that district, and minC is a VALID + MINIMAL length-s
    a,b-separator (removing minC leaves dist(a,b) > s).
  - COMPLETENESS: for every district that violates contiguity or the diameter
    bound, at least one cut with both endpoints in that district is produced
    (so the infeasible integer solution is always rejected). Conversely, a
    feasible solution produces no cuts.

This mirrors the branching of the original networkx callback and guarantees the
optimized separation accepts exactly the same set of integer solutions.
"""
import os
import json
import random
import networkx as nx
import fast_graph as fg


def dist_minus_set(G, a, b, removed):
    H = G.copy()
    H.remove_nodes_from([v for v in removed if v not in (a, b)])
    try:
        return nx.shortest_path_length(H, a, b)
    except nx.NetworkXNoPath:
        return float('inf')


def district_violates(G, Vj, s):
    """Reference: does district Vj violate contiguity or diameter<=s?"""
    if len(Vj) <= 1:
        return False
    H = G.subgraph(Vj)
    if not nx.is_connected(H):
        return True
    return nx.diameter(H) > s


def run_instance(G, k, s, selector, trials, seed):
    G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    rnd = random.Random(seed)

    minority_range = range(k)
    majority_range = range(0)
    if selector == 'both':
        kmin = max(1, k // 2)
        minority_range = range(kmin)
        majority_range = range(kmin, k)

    total_cuts = 0
    for _ in range(trials):
        # random assignment of every node to a district label in [0,k)
        assign = [rnd.randrange(k) for _ in range(n)]
        node_to_label = {v: assign[v] for v in range(n)}

        xval = {}
        yval = {}
        if selector == 'minority':
            for v in range(n):
                for j in range(k):
                    xval[v, j] = 1.0 if assign[v] == j else 0.0
        elif selector == 'majority':
            for v in range(n):
                for j in range(k):
                    yval[v, j] = 1.0 if assign[v] == j else 0.0
        else:  # both
            for v in range(n):
                for j in range(k):
                    if assign[v] == j and j in minority_range:
                        xval[v, j] = 1.0
                    else:
                        xval[v, j] = 0.0
                    if assign[v] == j and j in majority_range:
                        yval[v, j] = 1.0
                    else:
                        yval[v, j] = 0.0

        cuts = fg.separate_compactness(csr, s, selector, xval, yval,
                                       minority_range, majority_range, k)
        total_cuts += len(cuts)

        # district membership
        districts = {j: [v for v in range(n) if assign[v] == j] for j in range(k)}

        # ---- SOUNDNESS ----
        covered_labels = set()
        for (a, b, minC) in cuts:
            la, lb = node_to_label[a], node_to_label[b]
            assert la == lb, f"endpoints in different districts: {a},{b}"
            covered_labels.add(la)
            Vj = districts[la]
            Vset = set(Vj)
            assert all(c not in Vset for c in minC), "minC overlaps district"
            assert dist_minus_set(G, a, b, set(minC)) > s, \
                f"cut not a valid length-{s} separator: a={a} b={b} minC={minC}"
            for c in minC:
                assert dist_minus_set(G, a, b, set(minC) - {c}) <= s, \
                    "separator not minimal"

        # ---- COMPLETENESS ----
        for j in range(k):
            # only scanned labels matter for this selector
            if selector == 'minority' and j not in minority_range:
                continue
            if district_violates(G, districts[j], s):
                assert j in covered_labels, \
                    f"violating district {j} produced no cut (selector={selector})"
            # feasible districts may still legitimately produce no cut
        # if NO district violates, there must be zero cuts
        any_viol = any(district_violates(G, districts[j], s)
                       for j in range(k)
                       if not (selector == 'minority' and j not in minority_range))
        if not any_viol:
            assert len(cuts) == 0, "cuts produced for a feasible solution"
    return total_cuts


if __name__ == "__main__":
    random.seed(7)
    print("Random + grid graphs:")
    configs = []
    for i in range(12):
        nn = random.randint(10, 45)
        p = random.uniform(0.1, 0.4)
        configs.append((nx.gnp_random_graph(nn, p, seed=i), f"gnp#{i}"))
    for (r, c) in [(5, 5), (6, 8), (7, 9)]:
        configs.append((nx.grid_2d_graph(r, c), f"grid{r}x{c}"))

    for (G, label) in configs:
        Gi = nx.convert_node_labels_to_integers(G)
        if not nx.is_connected(Gi):
            continue
        for selector in ['minority', 'both']:
            for s in [2, 3, 4]:
                for k in [3, 5]:
                    tc = run_instance(Gi, k, s, selector, trials=6, seed=hash((label, s, k)) % 9999)
        print(f"  [{label}] OK")

    print("Real county instances:")
    base = "/sessions/gifted-bold-darwin/mnt/Sam_districting_paper/code/raw_data/county/json"
    for st in ["MS", "LA", "AL", "NM", "CO"]:
        path = os.path.join(base, f"{st}_counties.json")
        if not os.path.exists(path):
            continue
        d = json.load(open(path))
        G = nx.Graph()
        G.add_nodes_from(range(len(d["nodes"])))
        for u, nbrs in enumerate(d["adjacency"]):
            for e in nbrs:
                G.add_edge(u, e["id"])
        for selector in ['minority', 'both']:
            for s in [3, 5, 7]:
                for k in [4, 6]:
                    run_instance(G, k, s, selector, trials=8, seed=hash((st, s, k)) % 9999)
        print(f"  [{st}] OK")

    print("\nALL SEPARATION TESTS PASSED")
