"""
fetch_current_plan.py -- download a state's enacted congressional-district plan
and save it as raw_data/current_plan/<STATE>.shp, the file get_current_plan_stat()
expects (used only for the two "enacted plan" columns of the instance-info table).

Source: U.S. Census TIGER/Line 2022, 118th-Congress districts (cd118). These are
the post-2020-census ("2020 cycle") enacted maps, consistent with the 2020-census
tract/county geography used everywhere else in the paper. The files are clean
district polygons in EPSG:4269 (NAD83), matching the existing LA.shp / MS.shp.

Run on a machine with internet (needs geopandas):
    python fetch_current_plan.py AZ
    python fetch_current_plan.py AZ NV          (one or more states)
"""
import os
import sys

import geopandas as gpd

# USPS -> Census state FIPS
FIPS = {
    "AL": "01", "AK": "02", "AZ": "04", "AR": "05", "CA": "06", "CO": "08",
    "CT": "09", "DE": "10", "FL": "12", "GA": "13", "HI": "15", "ID": "16",
    "IL": "17", "IN": "18", "IA": "19", "KS": "20", "KY": "21", "LA": "22",
    "ME": "23", "MD": "24", "MA": "25", "MI": "26", "MN": "27", "MS": "28",
    "MO": "29", "MT": "30", "NE": "31", "NV": "32", "NH": "33", "NJ": "34",
    "NM": "35", "NY": "36", "NC": "37", "ND": "38", "OH": "39", "OK": "40",
    "OR": "41", "PA": "42", "RI": "44", "SC": "45", "SD": "46", "TN": "47",
    "TX": "48", "UT": "49", "VT": "50", "VA": "51", "WA": "53", "WV": "54",
    "WI": "55", "WY": "56",
}

# number of congressional (118th) districts per state -- used to verify that the
# downloaded enacted plan has the right number of districts
CD = {
    "AL": 7, "AK": 1, "AZ": 9, "AR": 4, "CA": 52, "CO": 8, "CT": 5, "DE": 1,
    "FL": 28, "GA": 14, "HI": 2, "ID": 2, "IL": 17, "IN": 9, "IA": 4, "KS": 4,
    "KY": 6, "LA": 6, "ME": 2, "MD": 8, "MA": 9, "MI": 13, "MN": 8, "MS": 4,
    "MO": 8, "MT": 2, "NE": 3, "NV": 4, "NH": 2, "NJ": 12, "NM": 3, "NY": 26,
    "NC": 14, "ND": 1, "OH": 15, "OK": 5, "OR": 6, "PA": 17, "RI": 2, "SC": 7,
    "SD": 1, "TN": 9, "TX": 38, "UT": 4, "VT": 1, "VA": 11, "WA": 10, "WV": 2,
    "WI": 8, "WY": 1,
}

OUT_DIR = os.path.join("raw_data", "current_plan")
URL = "https://www2.census.gov/geo/tiger/TIGER2022/CD/tl_2022_%s_cd118.zip"


def fetch(state):
    state = state.upper()
    if state not in FIPS:
        print("  unknown state:", state)
        return
    url = URL % FIPS[state]
    print("downloading", url)
    gdf = gpd.read_file(url)   # geopandas/pyogrio reads the zipped shapefile from the URL
    # Some cd118 files include a placeholder district coded 'ZZ' (water / areas not
    # assigned to a numbered district); drop it so only the real congressional
    # districts remain (e.g. IL: 18 -> 17).
    code_cols = [c for c in gdf.columns if c.upper().startswith("CD") and c.upper().endswith("FP")]
    for c in code_cols:
        gdf = gdf[gdf[c].astype(str).str.upper() != "ZZ"]
    # NAD83 (EPSG:4269), like the other states
    gdf = gdf.to_crs(epsg=4269)
    os.makedirs(OUT_DIR, exist_ok=True)
    dst = os.path.join(OUT_DIR, "%s.shp" % state)
    gdf.to_file(dst)
    print("  wrote %s  (%d districts, CRS %s)" % (dst, len(gdf), gdf.crs))
    return len(gdf)


def verify_all():
    """Re-download every state used in the paper (the 'paper' + 'hispanic_paper'
    instance sets) from the single authoritative cd118 source, and flag any whose
    district count does not match its congressional k."""
    import json
    with open("data.json") as f:
        data = json.load(f)
    states = []
    for ds in ("paper", "hispanic_paper"):
        for req in data.get(ds, []):
            s = req["state"].upper()
            if s not in states:
                states.append(s)
    print("Verifying %d paper states from the cd118 source: %s\n"
          % (len(states), ", ".join(states)))
    results = []
    for s in states:
        try:
            n = fetch(s)
        except Exception as e:
            print("  FAILED for %s: %s" % (s, e))
            n = None
        results.append((s, n, CD.get(s)))
    print("\n================ verification summary ================")
    bad = 0
    for s, n, k in results:
        if n is None:
            print("  %-3s  DOWNLOAD FAILED" % s); bad += 1
        elif k is None:
            print("  %-3s  %d districts  (no congressional k on record)" % (s, n))
        elif n == k:
            print("  %-3s  %d districts  == k  OK" % (s, n))
        else:
            print("  %-3s  %d districts  != k=%d  *** MISMATCH ***" % (s, n, k)); bad += 1
    print("======================================================")
    print("%d state(s) OK, %d flagged." % (len(results) - bad, bad))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage:")
        print("  python fetch_current_plan.py <STATE> [STATE ...]   e.g.  AZ IL")
        print("  python fetch_current_plan.py --verify-all          (re-download + check all paper states)")
        sys.exit(1)
    if sys.argv[1] == "--verify-all":
        verify_all()
    else:
        for s in sys.argv[1:]:
            try:
                fetch(s)
            except Exception as e:
                print("  FAILED for %s: %s" % (s, e))
                print("  (manual fallback: download %s , unzip, and rename the .shp/.shx/"
                      ".dbf/.prj to %s.* in raw_data/current_plan/)" % (URL % FIPS.get(s.upper(), "??"), s.upper()))
