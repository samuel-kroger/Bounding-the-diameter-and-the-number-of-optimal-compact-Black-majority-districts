@echo off
REM ===================================================================
REM  Table 4 -- Three-phase variable-fixing procedure
REM  Reports the diameter lower-bound improvement and the percentage of
REM  variables fixed per phase.  The radius-s ball is built with array BFS
REM  (no nx.ego_graph), and the combinatorial pre-screen avoids a MIP for
REM  most parcels.
REM  Output:  results_<group>\fixing_info.csv  (+ fixing_info_each_pass.csv)
REM ===================================================================
call conda activate base
cd /d "%~dp0"
python run_table.py --table fixing --group black    --dataset paper
python run_table.py --table fixing --group hispanic --dataset hispanic_paper
echo.
echo Wrote results_black\fixing_info.csv and results_hispanic\fixing_info.csv
pause
