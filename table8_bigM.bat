@echo off
REM ===================================================================
REM  Table 8 -- Separate-variable formulation vs. the big-M alternative
REM  (Reviewer 2, Minor comment 6).  Reports, for each instance and each
REM  formulation, the root LP-relaxation bound, the number of branch-and-
REM  bound nodes, and the solve time.  The separate-variable subpolytope is
REM  integral (Proposition 2), so its LP relaxation is tighter and it solves
REM  in fewer nodes / less time.
REM  Output:  bigM_comparison.csv
REM ===================================================================
call conda activate base
cd /d "%~dp0"
python run_bigM_comparison.py
echo.
echo Wrote bigM_comparison.csv
pause
