@echo off
REM ===================================================================
REM  Table 3 -- Upper bound on the number of majority-minority districts
REM  Solves the small upper-bounding MIP (Section 6); also reports the
REM  fraction of x/z/w variables fixed.  One-hour limit per instance.
REM  Output:  results_<group>\upper_bound_minority_table.csv
REM ===================================================================
call conda activate base
cd /d "%~dp0"
python run_table.py --table upper_bound_mm --group black    --dataset paper
python run_table.py --table upper_bound_mm --group hispanic --dataset hispanic_paper
echo.
echo Wrote results_black\upper_bound_minority_table.csv and results_hispanic\upper_bound_minority_table.csv
pause
