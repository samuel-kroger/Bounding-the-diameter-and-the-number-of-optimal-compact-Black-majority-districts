import os, json, random
import networkx as nx
import fast_graph as fg

def dist_minus(G,a,b,removed):
    H=G.copy(); H.remove_nodes_from([v for v in removed if v not in (a,b)])
    try: return nx.shortest_path_length(H,a,b)
    except nx.NetworkXNoPath: return float('inf')

def run(G,k,trials,seed):
    G=nx.convert_node_labels_to_integers(G); n=G.number_of_nodes()
    csr=fg.CSRGraph(n,list(G.edges())); rnd=random.Random(seed); minr=range(k)
    for _ in range(trials):
        assign=[rnd.randrange(k) for _ in range(n)]
        xval={(v,j):(1.0 if assign[v]==j else 0.0) for v in range(n) for j in range(k)}
        cuts=fg.separate_contiguity(csr,'minority',xval,{},minr,range(0),k)
        districts={j:[v for v in range(n) if assign[v]==j] for j in range(k)}
        covered=set()
        for a,b,C in cuts:
            assert assign[a]==assign[b]
            covered.add(assign[a])
            Vset=set(districts[assign[a]])
            assert all(c not in Vset for c in C)
            # C must fully separate a,b (contiguity = length infinity)
            assert dist_minus(G,a,b,set(C))==float('inf'), "not a separator"
        for j in range(k):
            Vj=districts[j]
            if len(Vj)>1 and not nx.is_connected(G.subgraph(Vj)):
                assert j in covered, f"disconnected district {j} not flagged"
        # feasible (all connected) -> no cuts
        if all(len(districts[j])<=1 or nx.is_connected(G.subgraph(districts[j])) for j in range(k)):
            assert len(cuts)==0
    return True

if __name__=="__main__":
    random.seed(3)
    for i in range(15):
        G=nx.gnp_random_graph(random.randint(8,40),random.uniform(0.1,0.4),seed=i)
        if not nx.is_connected(G): continue
        for k in [3,5]: run(G,k,6,seed=i*7+k)
    base="/sessions/gifted-bold-darwin/mnt/Sam_districting_paper/code/raw_data/county/json"
    for st in ["MS","AL","LA","CO","NM"]:
        d=json.load(open(os.path.join(base,f"{st}_counties.json")))
        G=nx.Graph(); G.add_nodes_from(range(len(d["nodes"])))
        for u,nbrs in enumerate(d["adjacency"]):
            for e in nbrs: G.add_edge(u,e["id"])
        for k in [4,6]: run(G,k,8,seed=hash(st)%999)
        print(f"  [{st}] OK")
    print("CONTIGUITY TESTS PASSED")
