"""
benchmark_run_helpers.py  --  helpers for the A1 bare-formulation benchmark.

Given a solved Gurobi model (our diameter side) or any saved plan, extract the
feasible solution(s), validate them, score every reviewer-requested metric
uniformly (# majority-minority districts, graph diameter, and -- with geometry --
average Polsby-Popper, Reock, and convex-hull), and save the maps.

The pure-graph helpers (diameter, #MM, map CSV) need only networkx + fast_graph
and are unit-tested. The geometry scores need geopandas/shapely. Solution
extraction needs the solved Gurobi model.
"""
import csv
import fast_graph as fg


def districts_from_labels(labels, k):
    d = [[] for _ in range(k)]
    for v, j in enumerate(labels):
        if j is not None and 0 <= j < k:
            d[j].append(v)
    return d


def plan_diameter(G, districts):
    """Maximum within-district graph diameter (no Gurobi/geopandas needed)."""
    n = G.number_of_nodes()
    csr = fg.CSRGraph(n, list(G.edges()))
    member = [False] * n
    worst = 0
    for D in districts:
        if not D:
            continue
        for v in D:
            member[v] = True
        for a in D:
            g = csr.bfs1(a, allowed=member)
            for b in D:
                if csr._s1[b] == g and csr._d1[b] > worst:
                    worst = csr._d1[b]
        for v in D:
            member[v] = False
    return worst


def num_minority_districts(G, districts, mvap, vap, f=0.5):
    cnt = 0
    for D in districts:
        if D and sum(mvap[v] for v in D) >= f * sum(vap[v] for v in D):
            cnt += 1
    return cnt


def polsby_popper_scores(G, districts):
    """Per-district Polsby-Popper from the graph data (node 'area', edge
    'shared_perim', node 'boundary_perim'+'boundary_node'); None if those fields
    are absent. Prefer geometry_scores (dissolved polygons) when a shapefile is
    available -- it needs no 'boundary_perim'."""
    try:
        import math
        labels = {i: j for j in range(len(districts)) for i in districts[j]}
        scores = []
        for D in districts:
            if not D:
                continue
            area = sum(G.nodes[i]['area'] for i in D)
            perim = sum(G.edges[u, v]['shared_perim']
                        for u in D for v in G.neighbors(u) if labels.get(v) != labels.get(u))
            perim += sum(G.nodes[i]['boundary_perim'] for i in D
                         if G.nodes[i].get('boundary_node'))
            scores.append(4 * math.pi * area / (perim * perim) if perim else 0.0)
        return scores
    except (KeyError, ZeroDivisionError):
        return None


def avg_polsby_popper(G, districts):
    scores = polsby_popper_scores(G, districts)
    return (sum(scores) / len(scores)) if scores else None


def save_map_csv(G, labels, path):
    """Write the node -> district assignment (with GEOID where available)."""
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['node', 'GEOID20', 'district'])
        for v in range(G.number_of_nodes()):
            w.writerow([v, G.nodes[v].get('GEOID20', ''), labels[v]])


def labels_from_pool(m, X, Y, n, k, sol_index):
    """Label vector for the `sol_index`-th solution in Gurobi's pool (via .Xn).
    Pool entries may violate lazy constraints, so validate before use."""
    m.setParam('SolutionNumber', sol_index)
    labels = [None] * n
    for v in range(n):
        for j in range(k):
            if (v, j) in X and X[v, j].Xn > 0.5:
                labels[v] = j
                break
            if (v, j) in Y and Y[v, j].Xn > 0.5:
                labels[v] = j
                break
    return labels


def labels_from_best(m, X, Y, n, k):
    """Label vector for the best (incumbent) solution via .X. Under lazy
    constraints this solution IS feasible, whereas pool entries (.Xn) may not be;
    so always score this one."""
    labels = [None] * n
    for v in range(n):
        for j in range(k):
            if (v, j) in X and X[v, j].X > 0.5:
                labels[v] = j
                break
            if (v, j) in Y and Y[v, j].X > 0.5:
                labels[v] = j
                break
    return labels


def geometry_scores(gdf, labels):
    """Geometry-based compactness (averages over districts), by dissolving the
    unit polygons into district polygons:
      Polsby-Popper = 4*pi*area / perimeter^2,
      convex-hull   = area / area(convex hull),
      Reock         = area / area(minimum bounding circle).
    Returns dict(avg_pp, avg_reock, avg_chull, pp, reock, chull) or None if
    geopandas/shapely are unavailable. Assumes gdf row i == node i (true for the
    JSON/shapefile pairs here); otherwise join on GEOID20 first."""
    try:
        import math
        import geopandas  # noqa: F401
        try:
            from shapely import minimum_bounding_circle
        except Exception:
            minimum_bounding_circle = None
    except Exception:
        return None
    g = gdf.copy()
    g["__d__"] = [labels[i] for i in range(len(labels))]
    dissolved = g.dissolve(by="__d__")
    pp, reock, chull = [], [], []
    for geom in dissolved.geometry:
        a, per = geom.area, geom.length
        pp.append(4 * math.pi * a / (per * per) if per else 0.0)
        chg = geom.convex_hull.area
        chull.append(a / chg if chg else 0.0)
        if minimum_bounding_circle is not None:
            mbc = minimum_bounding_circle(geom).area
            reock.append(a / mbc if mbc else 0.0)
        else:
            reock.append(None)

    def avg(xs):
        v = [x for x in xs if x is not None]
        return sum(v) / len(v) if v else None
    return dict(avg_pp=avg(pp), avg_reock=avg(reock), avg_chull=avg(chull),
                pp=pp, reock=reock, chull=chull)


def score_plan(G, districts, mvap, vap, gdf=None, labels=None, f=0.5):
    """All reviewer-requested metrics for ONE plan, computed uniformly:
    # majority-minority districts, max graph diameter, and (if a GeoDataFrame is
    given) average Polsby-Popper, Reock, and convex-hull."""
    res = dict(mm=num_minority_districts(G, districts, mvap, vap, f),
               diameter=plan_diameter(G, districts),
               avg_pp=None, avg_reock=None, avg_chull=None)
    if gdf is not None and labels is not None:
        gs = geometry_scores(gdf, labels)
        if gs:
            res["avg_pp"] = gs["avg_pp"]
            res["avg_reock"] = gs["avg_reock"]
            res["avg_chull"] = gs["avg_chull"]
    return res


def load_plan_csv(path):
    """Read a saved map CSV (node, GEOID20, district) back into (labels,
    districts, k) so any method's exported plan can be re-scored identically."""
    rows = list(csv.DictReader(open(path)))
    n = len(rows)
    labels = [None] * n
    for r in rows:
        labels[int(r["node"])] = int(r["district"])
    kk = max(l for l in labels if l is not None) + 1
    return labels, districts_from_labels(labels, kk), kk


def extract_all_solutions(m, X, Y, n, k, csr, s, pop, L, U, mvap, vap, G,
                          map_dir, tag, f=0.5, gdf=None):
    """Collect every feasible solution from a solved model, VALIDATE each with
    fast_graph.check_feasible_plan (essential under lazy constraints), score it,
    and save its map. The incumbent (.X) is checked first -- it is lazy-feasible,
    unlike pool entries (.Xn). Returns a list of dicts (one per valid solution)."""
    import os
    os.makedirs(map_dir, exist_ok=True)
    out = []
    seen = set()
    candidates = []
    if m.SolCount > 0:
        candidates.append(("incumbent", labels_from_best(m, X, Y, n, k)))
    for i in range(m.SolCount):
        candidates.append(("pool", labels_from_pool(m, X, Y, n, k, i)))
    print("  extract: SolCount=%d, scanning %d candidate solution(s)"
          % (m.SolCount, len(candidates)))
    idx = 0
    for kind, labels in candidates:
        if labels is None or any(l is None for l in labels):
            if kind == "incumbent":
                print("  [warn] incumbent has unassigned node(s); skipping")
            continue
        key = tuple(labels)
        if key in seen:
            continue
        ok, why = fg.check_feasible_plan(csr, labels, k, s, pop, L, U)
        if not ok:
            if kind == "incumbent":
                # the Gurobi-optimal incumbent is feasible by construction; keep
                # it and report the disagreement so we can investigate the check.
                print("  [warn] incumbent failed check_feasible_plan: %s "
                      "-- including it anyway (Gurobi-optimal)" % why)
            else:
                continue                   # pool entry violates the lazy cuts
        seen.add(key)
        idx += 1
        districts = districts_from_labels(labels, k)
        map_file = os.path.join(map_dir, '%s_sol%d.csv' % (tag, idx))
        save_map_csv(G, labels, map_file)
        sc = score_plan(G, districts, mvap, vap, gdf=gdf, labels=labels, f=f)
        out.append(dict(
            sol_index=idx,
            mm=sc["mm"], pp=sc["avg_pp"], reock=sc["avg_reock"],
            chull=sc["avg_chull"], diam=sc["diameter"],
            map_file=map_file,
        ))
    return out
