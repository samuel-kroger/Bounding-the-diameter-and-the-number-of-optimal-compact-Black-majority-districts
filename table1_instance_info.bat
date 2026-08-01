@echo off
REM ===================================================================
REM  Table 1 -- Instance information
REM  Columns: state, level, k, n, m, graph diameter, % minority VAP, ...
REM  Output:  results_black\instance_info.csv  and  results_hispanic\instance_info.csv
REM ===================================================================
call conda activate base
cd /d "%~dp0"
python run_table.py --table instance_info --group black    --dataset paper
python run_table.py --table instance_info --group hispanic --dataset hispanic_paper
echo.
echo Wrote results_black\instance_info.csv and results_hispanic\instance_info.csv
pause
