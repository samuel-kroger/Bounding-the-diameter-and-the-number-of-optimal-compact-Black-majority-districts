@echo off
REM ===================================================================
REM  Table 5, Panel B -- Symmetry-breaking in the FEASIBLE regime
REM
REM  Panel A (table5_symmetry.bat) runs at s = ell_s, where the tract
REM  instances have NO valid districting: it measures how fast infeasibility
REM  can be PROVEN with vs without the symmetry-breaking constraints.
REM
REM  Panel B (this file) runs at an s where a plan is KNOWN to exist (taken
REM  from the A1 benchmark's feasible solutions), so both arms have a feasible
REM  region. That is what populates obj / bound / gap / gap_improvement.
REM
REM  Parcel fixing is ON in both arms (as in the original symmetry experiment),
REM  so the only difference between the two columns is the symmetry constraints.
REM  No warm start is used -- a warm start disables symmetry by design.
REM
REM  Instances (s = min max-diameter among our feasible solutions):
REM    MS county s=6 (#MM=1)   MS tract s=16 (#MM=1)
REM    SC tract  s=16 (#MM=1)  LA tract  s=33 (#MM=1)
REM    NM tract  s=16 (#MM=2, hispanic)
REM
REM  Output: results_<group>\symmetry_table_feasible.csv
REM  Resumable: rows are keyed on (state, parcel level, s), so re-running
REM  after an interruption skips what is already done.
REM ===================================================================
call conda activate base
cd /d "%~dp0"
python run_table.py --table symmetry_feasible --group black    --dataset symmetry_feasible
python run_table.py --table symmetry_feasible --group hispanic --dataset symmetry_feasible_hispanic
echo.
echo Wrote results_black\symmetry_table_feasible.csv and results_hispanic\symmetry_table_feasible.csv
pause
