"""
repair_connectivity.py -- reconnect the FL and NY tract contiguity graphs.

WHY
---
Three of the states with the largest Latino populations (CA, FL, NY) are absent
from the paper's tract-level experiments because their RAW tract contiguity
graphs are disconnected. This is a property of the census geography (islands and
offshore water tracts), not of the formulation. Following Validi et al. (2022),
who add contiguity edges across the natural geographic discontinuities (bridges
and ferry routes), we repair the graphs by adding a small number of explicit
water-crossing edges.

The disconnection is tiny. Each of FL and NY has exactly ONE stranded tract:

  FL  node 1636  GEOID 12087980100  Monroe County tract 9801   pop     11
      The Dry Tortugas: a national-park tract about 110 km west of Key West,
      reachable only by boat or seaplane. Its nearest neighbor in the graph is
      the Monroe County offshore water tract 9900 (node 51, pop 0), which is
      already a degree-25 hub joining the Keys. We also add the edge to the
      nearest POPULATED Monroe tract, 9725 (node 2077), which is the Key West
      end of the Dry Tortugas ferry.

  NY  node 2678  GEOID 36061000100  New York County tract 1    pop      0
      Governors Island. Ferries run to it from Lower Manhattan (Battery
      Maritime Building) and from Red Hook in Brooklyn, so we add one edge for
      each: Manhattan tract 5 (node 2681) and Kings County tract 53.01
      (node 3267). Both are ~1.05 km away, the two nearest tracts in the state.

CA is repaired too. It is the largest instance in the study (n = 9,129, k = 52,
against TX at n = 6,896, k = 38), so expect it to be the slowest by some margin.

WHAT IT DOES
------------
Rewrites raw_data/tract/json/{FL,NY}_tracts.json with the extra edges, keeping a
one-time backup at *_tracts.ORIGINAL.json. It is idempotent: re-running detects
that the edges are already present and does nothing. Node attributes, node ids,
and node ordering are preserved exactly, so every downstream table is unaffected
apart from the added adjacencies.

USAGE
-----
    python repair_connectivity.py            # repair FL and NY
    python repair_connectivity.py --check    # report status, change nothing
    python repair_connectivity.py --restore  # put the original graphs back
    python repair_connectivity.py --auto     # recompute nearest tracts instead
                                             # of using the reviewed edge list
"""
import argparse
import json
import os
import sys

import networkx as nx
from networkx.readwrite import json_graph

TRACT_JSON = "./raw_data/tract/json/{st}_tracts.json"
TRACT_SHP = "./raw_data/tract/shape/{st}_tracts.shp"
BACKUP = "./raw_data/tract/json/{st}_tracts.ORIGINAL.json"

# Reviewed water-crossing edges, keyed by state. Each entry is
#   (node_u, node_v, "justification")
# Node ids are positional and match both the JSON and the shapefile row order.
REPAIRS = {
    "FL": [
        (1636, 51, "Dry Tortugas (tract 9801) to Monroe offshore water tract 9900"),
        (1636, 2077, "Dry Tortugas to Key West tract 9725, the ferry origin"),
    ],
    "NY": [
        (2678, 2681, "Governors Island to Manhattan tract 5, Battery ferry"),
        (2678, 3267, "Governors Island to Brooklyn tract 53.01, Red Hook ferry"),
    ],
    # CA has three stranded tracts in two components, all of them islands.
    #   1212  06075980401  Farallon Islands, San Francisco County, pop 0.
    #         Uninhabited; administratively part of San Francisco, so we attach
    #         it to San Francisco's own offshore water tract 9901 (node 2346).
    #   3181  06037599000  Avalon, Santa Catalina Island, LA County, pop 3,322.
    #   7374  06037599100  the rest of Santa Catalina plus the other LA County
    #         Channel Islands, pop 553. These two are already adjacent to each
    #         other, forming the size-2 component, so one link to the mainland
    #         suffices; we add both the water-tract link and the ferry link.
    #         The Catalina Express sails from San Pedro, tract 2975.01.
    "CA": [
        (1212, 2346, "Farallon Islands to San Francisco water tract 9901"),
        (3181, 2441, "Santa Catalina (Avalon) to LA County water tract 9903"),
        (3181, 3840, "Santa Catalina (Avalon) to San Pedro tract 2975.01, ferry"),
        (7374, 2441, "LA County Channel Islands to LA water tract 9903"),
    ],
}

# Sanity guards: the repair refuses to run if the geography moved under it.
EXPECTED_GEOID = {
    "FL": {1636: "12087980100", 51: "12087990000", 2077: "12087972500"},
    "NY": {2678: "36061000100", 2681: "36061000500", 3267: "36047005301"},
    "CA": {1212: "06075980401", 2346: "06075990100", 3181: "06037599000",
           2441: "06037990300", 3840: "06037297501", 7374: "06037599100"},
}


def load(st):
    with open(TRACT_JSON.format(st=st)) as fh:
        data = json.load(fh)
    return data, json_graph.adjacency_graph(data, multigraph=False)


def describe(G, v):
    a = G.nodes[v]
    return "node {:5d}  GEOID {}  tract {:>8}  pop {:>6}".format(
        v, a.get("GEOID20", "?"), a.get("NAME20", "?"), a.get("P0010001", "?"))


def components(G):
    return sorted((sorted(c) for c in nx.connected_components(G)), key=len, reverse=True)


def check(st):
    _, G = load(st)
    comps = components(G)
    ok = len(comps) == 1
    print("{}  n = {:,}  m = {:,}  components = {}  {}".format(
        st, G.number_of_nodes(), G.number_of_edges(), len(comps),
        "CONNECTED" if ok else "DISCONNECTED"))
    for c in comps[1:]:
        for v in c:
            print("      stranded:", describe(G, v))
    return ok


def nearest_in_main(st, G, v, how_many=2):
    """Geometric fallback: the nearest tracts in the main component, by polygon
    distance in Web Mercator. Used only with --auto."""
    import geopandas as gpd
    gdf = gpd.read_file(TRACT_SHP.format(st=st))
    if gdf.crs is None:
        gdf = gdf.set_crs(4269)
    gdf = gdf.to_crs(3857)
    main = set(components(G)[0])
    dist = gdf.geometry.distance(gdf.geometry.iloc[v])
    order = [i for i in dist.sort_values().index if i in main]
    return order[:how_many]


def verify_geoids(st, G):
    for node, geoid in EXPECTED_GEOID[st].items():
        actual = G.nodes[node].get("GEOID20")
        if actual != geoid:
            sys.exit(
                "ABORT: {} node {} is GEOID {}, expected {}. The raw data has "
                "changed; re-derive the repair edges before running."
                .format(st, node, actual, geoid))


def repair(st, auto=False):
    path = TRACT_JSON.format(st=st)
    backup = BACKUP.format(st=st)
    data, G = load(st)

    comps = components(G)
    if len(comps) == 1:
        print("{}: already connected, nothing to do.".format(st))
        return False

    print("{}: {} components, repairing.".format(st, len(comps)))

    if auto:
        edges = []
        for c in comps[1:]:
            for v in c:
                for u in nearest_in_main(st, G, v, 2):
                    edges.append((v, u, "auto: nearest tract by polygon distance"))
    else:
        verify_geoids(st, G)
        edges = REPAIRS[st]

    added = 0
    for u, v, why in edges:
        if G.has_edge(u, v):
            print("      already present: {} -- {}".format(u, v))
            continue
        G.add_edge(u, v)
        added += 1
        print("      + edge {:5d} -- {:<5d}  {}".format(u, v, why))
        print("            {}".format(describe(G, u)))
        print("            {}".format(describe(G, v)))

    if not nx.is_connected(G):
        sys.exit("ABORT: {} still disconnected after repair.".format(st))

    if not os.path.exists(backup):
        with open(backup, "w") as fh:
            json.dump(data, fh)
        print("      backup written: {}".format(backup))

    out = json_graph.adjacency_data(G)
    # keep the original top-level flags so the file stays byte-compatible with
    # whatever loader is used downstream (gerrychain or plain networkx)
    for key in ("directed", "multigraph", "graph"):
        if key in data:
            out[key] = data[key]
    with open(path, "w") as fh:
        json.dump(out, fh)

    print("      {} edges added, now connected. m: {:,} -> {:,}".format(
        added, G.number_of_edges() - added, G.number_of_edges()))
    return True


def restore(st):
    backup = BACKUP.format(st=st)
    if not os.path.exists(backup):
        print("{}: no backup to restore.".format(st))
        return
    with open(backup) as fh:
        data = json.load(fh)
    with open(TRACT_JSON.format(st=st), "w") as fh:
        json.dump(data, fh)
    print("{}: original graph restored.".format(st))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--states", nargs="+", default=["FL", "NY", "CA"])
    ap.add_argument("--check", action="store_true", help="report only")
    ap.add_argument("--restore", action="store_true", help="undo the repair")
    ap.add_argument("--auto", action="store_true",
                    help="recompute nearest tracts instead of the reviewed list")
    args = ap.parse_args()

    if args.restore:
        for st in args.states:
            restore(st)
        return

    if args.check:
        allok = all(check(st) for st in args.states)
        sys.exit(0 if allok else 1)

    for st in args.states:
        repair(st, auto=args.auto)
    print()
    for st in args.states:
        check(st)


if __name__ == "__main__":
    main()
