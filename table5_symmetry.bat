@echo off
REM ===================================================================
REM  Table 5 -- Symmetry-breaking speedup
REM  Solves each instance with and without the independent-set symmetry
REM  constraints and reports the speedup (e.g., ~200x on MS county).
REM  Output:  results_<group>\symmetry_table.csv
REM
REM  Run AFTER table2 (lower_bound_s) and table4 (fixing): this reuses the
REM  cached ell_s and the independent-set / fixing files.
REM ===================================================================
call conda activate base
cd /d "%~dp0"
python run_table.py --table symmetry --group black    --dataset paper
python run_table.py --table symmetry --group hispanic --dataset hispanic_paper
echo.
echo Wrote results_black\symmetry_table.csv and results_hispanic\symmetry_table.csv
pause
