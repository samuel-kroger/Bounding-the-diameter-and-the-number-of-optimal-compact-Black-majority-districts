@echo off
REM ===================================================================
REM  Table 7 -- Benchmark: our diameter model vs. Belotti (PP-optimal MISOCP
REM  and Carve-Complete-Cleanup) and GerryChain, on a common instance set.
REM
REM  PREREQUISITES (generate the competitors' plans once; see RUN_ALL_README.md):
REM     cd ..\Polsby_Popper_optimization-main
REM     python build_belotti_data.py            (Belotti-format data from our graphs)
REM     run_belotti_our_instances.bat           (Belotti PP-optimal, county)
REM     run_belotti_tracts.bat                  (Belotti PP-optimal, tracts)
REM     run_gingles.bat                         (Belotti CCC, MM heuristic)
REM     cd ..\code  &  python run_gerrychain.py (GerryChain upper bounds + warm starts)
REM
REM  This batch (re)runs OUR exact model with the FINAL code (warm-started from
REM  the GerryChain/CCC plans) and rebuilds the combined comparison workbook.
REM  Output:  A1_benchmark_results.xlsx  and  A1_benchmark_combined.xlsx
REM ===================================================================
call conda activate base
cd /d "%~dp0"
echo [1/2] Our exact diameter model (warm-started)...
python run_A1_ours.py
echo [2/2] Scoring all methods into the combined workbook...
python score_belotti_plans.py
echo.
echo Wrote A1_benchmark_results.xlsx and A1_benchmark_combined.xlsx
pause
