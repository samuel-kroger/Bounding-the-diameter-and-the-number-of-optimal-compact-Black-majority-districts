"""
make_tradeoff_figure.py -- redraw Figure 1 of the paper (tradeoff_frontier.pdf).

The repository had no script for this figure, so this one reproduces the
published style exactly and regenerates it from results_tradeoff/
tradeoff_frontier.csv, which now includes the FL, NY and CA rows.

STYLE (matching the published figure)
  * two panels, (a) Majority-Black and (b) Majority-Latino
  * instances that reach 2 or more majority-minority districts are drawn as a
    staircase, x = diameter bound s, y = number of districts, steps-post, with
    a marker at each level and the state named in the panel legend
  * instances that reach exactly 1 are a single gray cross: one point is not a
    curve, and naming a dozen of them would swamp the legend
  * the two county-level instances whose minimum diameter is CERTIFIED by the
    exact model are red stars, annotated
  * instances with no positive majority-minority incumbent are omitted, and
    the count is reported so the caption can state it. Note this is "none
    found", not "none exists": several have a positive k_m.

USAGE
    python make_tradeoff_figure.py                 # writes ../paper/tradeoff_frontier.pdf
    python make_tradeoff_figure.py --out fig.pdf
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

CSV = "results_tradeoff/tradeoff_frontier.csv"
DEFAULT_OUT = "../paper/tradeoff_frontier.pdf"

PANELS = [("black", "(a) Majority-Black"), ("hispanic", "(b) Majority-Latino")]

# One colour per state, held FIXED across both panels. Several states appear in
# both (NY, IL, TX, NJ), and letting matplotlib's cycler assign per-panel makes
# the same state a different colour on the left and the right, which invites
# misreading. A fixed key removes that.
# Within a panel every colour must be distinct; across panels a state keeps its
# colour. CA gets the strong red because its curve is the headline result, so NY
# (which appears in BOTH panels) takes brown rather than colliding with it.
STATE_COLOR = {
    # majority-Black panel
    "GA": "#1f77b4", "MD": "#2ca02c", "IL": "#ff7f0e", "NY": "#8c564b",
    # majority-Latino panel
    "AZ": "#1f77b4", "NM tract": "#ff7f0e", "TX": "#2ca02c",
    "CA": "#d62728", "FL": "#9467bd",
}

# Some frontiers coincide exactly. TX reaches all ten of its levels at s = 20 and
# CA reaches eleven of its fourteen at the same s, so the two curves lie on top
# of one another and whichever is drawn second hides the first completely.
# Dashing the shorter one lets both read at the same x.
STATE_DASH = {"TX": (4, 2.4), "AZ": (6, 2)}

# A few labels would run into a neighbouring curve if placed to the right of
# their end point, so they hang to the left instead.
LABEL_LEFT = {"NM tract"}

# Per-label nudges, in typographic points, keyed by (panel, label). The default
# placement sits a curve's name up and to the right of its top point, which
# collides where curves crowd together. These overrides move individual labels
# clear of neighbouring lines and markers.
#   dx > 0 places the label to the RIGHT of the point, dx < 0 to the LEFT
#   dy > 0 raises it, dy < 0 lowers it
LABEL_OFFSET = {
    ("black", "NY"): (-4, 4),      # left, clear of the IL and MD verticals
    ("black", "IL"): (-11, -7),    # left and down, clear of the NY vertical at s=21
    ("black", "MD"): (7, -9),      # down, away from the GA label above it
    ("black", "GA"): (7, -9),      # down, tucked beside its own top marker
    ("hispanic", "CA"): (7, -9),   # down, off the top frame
    ("hispanic", "AZ"): (-6, 2),   # left of the vertical, not on top of it
    ("hispanic", "NY"): (6, -11),  # down and just right: the space to the LEFT is
                                   # occupied by the CA/TX column at s = 20, so a
                                   # leftward label would sit on top of it
}


def parse_frontier(text):
    """'1:19 2:19 3:19 4:24' -> ([19,19,19,24], [1,2,3,4])"""
    xs, ys = [], []
    for tok in (text or "").split():
        if ":" not in tok:
            continue
        k, s = tok.split(":")
        ys.append(int(k))
        xs.append(int(s))
    pairs = sorted(zip(ys, xs))
    return [p[1] for p in pairs], [p[0] for p in pairs]


def label_for(row):
    """County-level instances share a state name with their tract instance, so
    disambiguate exactly as the published figure does ('NM tract')."""
    st, lvl = row["State"].strip(), row["level"].strip()
    return st if lvl == "tract" and st not in ("NM", "MS") else "%s %s" % (st, lvl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--csv", default=CSV)
    args = ap.parse_args()

    with open(args.csv) as fh:
        rows = [r for r in csv.DictReader(fh) if r.get("State")]

    fig, axes = plt.subplots(1, 2, figsize=(10.88, 4.55))
    omitted = {"black": [], "hispanic": []}

    for ax, (grp, title) in zip(axes, PANELS):
        mine = [r for r in rows if r["group"].strip() == grp]
        singles_x = []
        curves = []
        curve_specs = []
        n_curves = 0

        for r in mine:
            xs, ys = parse_frontier(r["frontier s_min(k_m)"])
            best = int(r["max_MM_found"] or 0)
            certified = (r.get("s*_proven") or "").strip().lower() == "yes"

            if best == 0 or not xs:
                omitted[grp].append(label_for(r))
                continue

            if certified:
                ax.plot(xs[0], ys[0], marker="*", markersize=17,
                        color="crimson", markeredgecolor="black",
                        linestyle="none", zorder=5)
                ax.annotate("%s ($s^*$)" % label_for(r),
                            (xs[0], ys[0]), textcoords="offset points",
                            xytext=(6, 10), fontsize=9)
                continue

            if best == 1:
                singles_x.append(xs[0])
                continue

            lab = label_for(r)
            curve_specs.append((max(ys), xs, ys, lab))
            # Label each curve at its top point rather than relying on a legend.
            # Curves here are mostly near-vertical and overlap heavily, so a
            # legend forces the reader to match colours; a label sitting on the
            # curve does not.
            n_curves += 1

        # Draw the TALLEST curve first so shorter ones land on top and stay
        # visible where they coincide.
        for _h, xs, ys, lab in sorted(curve_specs, key=lambda c: -c[0]):
            col = STATE_COLOR.get(lab)
            kw = {}
            if lab in STATE_DASH:
                kw["dashes"] = STATE_DASH[lab]
            ax.plot(xs, ys, drawstyle="steps-post", marker="o",
                    markersize=4, linewidth=1.9, color=col, zorder=3, **kw)
            curves.append((xs[-1], ys[-1], lab, col))

        if singles_x:
            ax.plot(singles_x, [1] * len(singles_x), marker="x",
                    linestyle="none", color="0.55", markersize=9,
                    markeredgewidth=1.6, zorder=2)

        ax.set_title(title, fontsize=12)
        ax.set_xlabel("diameter bound $s$")
        ax.set_ylabel("number of majority-minority districts")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        # district counts are integers; without this matplotlib labels the
        # short axis in half-districts (0.5, 1.5, ...)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        # headroom so the end labels are not clipped at the frame
        ax.margins(x=0.11, y=0.13)
        # place the end labels, nudging apart any that would collide
        curves.sort(key=lambda c: (-c[1], c[0]))
        used = []
        for x, y, lab, col in curves:
            override = LABEL_OFFSET.get((grp, lab))
            if override is not None:
                dx, dy = override
                left = dx < 0
            else:
                left = lab in LABEL_LEFT
                dx, dy = (-7, 4) if left else (7, 4)
                # only auto-nudge labels that have no explicit placement
                for ux, uy in used:
                    if abs(ux - x) < 1.6 and abs(uy - y) < 1.2:
                        dy += 13
            ax.annotate(lab, (x, y), textcoords="offset points",
                        xytext=(dx, dy), fontsize=9.5, color=col,
                        fontweight="bold", zorder=6,
                        ha="right" if left else "left")
            used.append((x, y))

    handles = [
        Line2D([], [], marker="*", color="crimson", markeredgecolor="black",
               markersize=15, linestyle="none",
               label="certified minimum $s^*$ (exact)"),
        Line2D([], [], marker="x", color="0.55", markersize=9,
               markeredgewidth=1.6, linestyle="none",
               label="upper bound (heuristic)"),
        Line2D([], [], marker="o", color="black", markersize=4, linewidth=1.6,
               label="multi-district frontier (heuristic)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=10)

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(args.out, bbox_inches="tight")
    print("wrote %s" % os.path.abspath(args.out))

    total = 0
    for grp, names in omitted.items():
        if names:
            print("  omitted (%s, no positive incumbent found): %s"
                  % (grp, ", ".join(sorted(names))))
            total += len(names)
    print("  TOTAL OMITTED = %d  <- the caption must state this number" % total)


if __name__ == "__main__":
    main()
