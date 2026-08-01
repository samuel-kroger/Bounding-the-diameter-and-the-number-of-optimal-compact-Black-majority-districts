"""
benchmark_to_excel.py  --  write the A1 benchmark results to a formatted Excel
workbook with two sheets:

  1. "A1 benchmark"            -- one row per (instance, method): the best result
                                  of each method (# MM districts, time, gap,
                                  avg Polsby-Popper, max diameter). Methods:
                                    * diameter (ours, bare)
                                    * PP-MISOCP (theirs, bare)        [min avg inv-PP, fixed k]
                                    * PP-MISOCP (theirs, referee variant)  [max #MM s.t. PP>=tau]
  2. "All feasible solutions"  -- one row per feasible solution found within the
                                  1-hour time limit, with its objective, time
                                  found, gap, PP, diameter, and saved map path.

All formulations are BARE (no fixing / symmetry / clique cuts / warm start /
heuristics) and solved with identical Gurobi settings (TimeLimit = 3600 s).

Library use:
    from benchmark_to_excel import build_workbook
    build_workbook(summary_rows, solution_rows, "A1_benchmark_results.xlsx")
Run with no args to emit a ready-to-fill template.
"""
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

FONT = "Arial"
NAVY = "1F3864"; BLUE = "2E5496"; WHITE = "FFFFFF"
SHADE_A = "FFFFFF"; SHADE_B = "EDF2FB"   # alternating per-instance group shades
METHOD_FILL = {"diameter (ours, full)": "E2EFDA",
               "Belotti et al. (published, PP-optimal)": "FCE4D6"}
_thin = Side(style="thin", color="BFBFBF")
BORDER = Border(left=_thin, right=_thin, top=_thin, bottom=_thin)

COLS = [
    ("state", "State", "@"), ("level", "Level", "@"),
    ("n", "n", "#,##0"), ("k", "k", "0"),
    ("method", "Method", "@"),
    ("compact", "compactness\n(s* or PP≥τ)", "@"),
    ("mm", "# MM", "0"), ("time", "time (s)", "#,##0.0"),
    ("gap", "gap", "0.0%"), ("pp", "avg PP", "0.000"),
    ("reock", "avg Reock", "0.000"), ("chull", "avg CH", "0.000"),
    ("diam", "max diam", "0"),
]
SOL_COLS = [("state", "State", "@"), ("level", "Level", "@"), ("method", "Method", "@"),
            ("sol_index", "Sol #", "0"), ("mm", "# MM", "0"),
            ("time_found", "time found (s)", "#,##0.0"), ("gap", "gap at end", "0.0%"),
            ("pp", "avg PP", "0.000"), ("reock", "avg Reock", "0.000"),
            ("chull", "avg CH", "0.000"), ("diam", "max diam", "0"),
            ("map_file", "map file", "@")]


def _summary_sheet(ws, rows):
    ws.title = "A1 benchmark"
    ncol = len(COLS)
    title = ("A1 benchmark — our diameter model (full, all accelerations) vs. "
             "Belotti, Buchanan & Ezazipour's published Polsby–Popper results")
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=ncol)
    c = ws.cell(1, 1, title); c.font = Font(name=FONT, bold=True, size=12, color=WHITE)
    c.fill = PatternFill("solid", fgColor=NAVY)
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.row_dimensions[1].height = 30
    note = ("Ours: the full optimized model (reduced formulation + parcel fixing + "
            "independent-set symmetry-breaking + clique distance cuts), maximizing #MM "
            "subject to diameter ≤ s, starting at s = ℓ_s and increasing until feasible "
            "(s*), TimeLimit = 3600 s. Theirs: Belotti et al.'s published Polsby-Popper "
            "results (their method optimizes PP) — filled from their paper/repo; their "
            "code is not re-run here.")
    ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=ncol)
    c = ws.cell(2, 1, note); c.font = Font(name=FONT, italic=True, size=9, color="404040")
    c.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    ws.row_dimensions[2].height = 44
    for j, (key, header, fmt) in enumerate(COLS, start=1):
        cell = ws.cell(3, j, header); cell.font = Font(name=FONT, bold=True, size=10, color=WHITE)
        cell.fill = PatternFill("solid", fgColor=BLUE)
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BORDER
    ws.row_dimensions[3].height = 30
    # group rows by instance (state, level) for alternating shade
    group_idx, prev = -1, None
    for i, row in enumerate(rows):
        key = (row.get("state"), row.get("level"))
        if key != prev:
            group_idx += 1; prev = key
        shade = SHADE_A if group_idx % 2 == 0 else SHADE_B
        r = 4 + i
        for j, (k_, header, fmt) in enumerate(COLS, start=1):
            cell = ws.cell(r, j, row.get(k_)); cell.font = Font(name=FONT, size=10)
            cell.number_format = fmt
            cell.alignment = Alignment(horizontal="left" if j <= 2 or j == 5 else "center",
                                       vertical="center")
            cell.border = BORDER
            if k_ == "method" and row.get("method") in METHOD_FILL:
                cell.fill = PatternFill("solid", fgColor=METHOD_FILL[row["method"]])
            else:
                cell.fill = PatternFill("solid", fgColor=shade)
    ws.freeze_panes = "A4"
    for j, w in enumerate([7, 8, 7, 5, 34, 13, 7, 9, 7, 8, 8, 8, 9], start=1):
        ws.column_dimensions[get_column_letter(j)].width = w


def _solutions_sheet(ws, solutions):
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(SOL_COLS))
    c = ws.cell(1, 1, "All feasible solutions found within the 1-hour time limit, with its map")
    c.font = Font(name=FONT, bold=True, size=12, color=WHITE); c.fill = PatternFill("solid", fgColor=NAVY)
    c.alignment = Alignment(horizontal="center", vertical="center"); ws.row_dimensions[1].height = 24
    for j, (k_, h, f) in enumerate(SOL_COLS, start=1):
        cell = ws.cell(2, j, h); cell.font = Font(name=FONT, bold=True, size=10, color=WHITE)
        cell.fill = PatternFill("solid", fgColor=BLUE)
        cell.alignment = Alignment(horizontal="center", vertical="center"); cell.border = BORDER
    for i, s in enumerate(solutions):
        r = 3 + i
        for j, (k_, h, f) in enumerate(SOL_COLS, start=1):
            cell = ws.cell(r, j, s.get(k_)); cell.font = Font(name=FONT, size=10); cell.number_format = f
            cell.alignment = Alignment(horizontal="center" if j > 3 else "left", vertical="center"); cell.border = BORDER
            if i % 2 == 1:
                cell.fill = PatternFill("solid", fgColor="F2F2F2")
    ws.freeze_panes = "A3"
    for j, w in enumerate([7, 8, 34, 6, 7, 13, 11, 8, 8, 8, 9, 34], start=1):
        ws.column_dimensions[get_column_letter(j)].width = w


def build_workbook(summary_rows, solution_rows, path="A1_benchmark_results.xlsx"):
    wb = Workbook()
    _summary_sheet(wb.active, summary_rows)
    _solutions_sheet(wb.create_sheet("All feasible solutions"), solution_rows)
    wb.save(path)
    return path


METHODS = ["diameter (ours, full)",
           "Belotti et al. (published, PP-optimal)"]
_INSTANCES = [("MS", "county", 82, 4), ("AL", "county", 67, 7),
              ("LA", "county", 64, 6), ("SC", "county", 46, 7),
              ("NM", "tract", 612, 3), ("MS", "tract", 878, 4),
              ("SC", "tract", 1323, 7), ("LA", "tract", 1388, 6),
              ("CO", "tract", 1447, 8)]
TEMPLATE_ROWS = [dict(state=s, level=l, n=n, k=k, method=m)
                 for (s, l, n, k) in _INSTANCES for m in METHODS]
TEMPLATE_SOLUTIONS = [
    dict(state="MS", level="county", method="diameter (ours, bare)", sol_index=1,
         map_file="maps/MS_county_ours_sol1.csv"),
]

if __name__ == "__main__":
    build_workbook(TEMPLATE_ROWS, TEMPLATE_SOLUTIONS, "A1_benchmark_results.xlsx")
    print("wrote A1_benchmark_results.xlsx (long-format summary with 3 methods + solutions)")
