@echo off
REM ===================================================================
REM  Table 2 -- Lower bound on the diameter parameter s (ell_s)
REM  Computed via the independence-number / power-graph certificate
REM  (array-based, no networkx), then binary search.
REM  Output:  results_<group>\lower_bound_s_table.csv
REM
REM  Run this BEFORE the upper-bound, fixing, and symmetry tables: the
REM  instance constructor caches ell_s, which those tables reuse.
REM ===================================================================
call conda activate base
cd /d "%~dp0"
python run_table.py --table lower_bound_s --group black    --dataset paper
python run_table.py --table lower_bound_s --group hispanic --dataset hispanic_paper
echo.
echo Wrote results_black\lower_bound_s_table.csv and results_hispanic\lower_bound_s_table.csv
pause
