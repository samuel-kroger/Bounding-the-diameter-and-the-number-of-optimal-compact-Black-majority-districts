@echo off
REM ============================================================
REM  A1 bare-formulation benchmark -- OUR diameter model side.
REM  Open "Anaconda Prompt", then run:   run_A1.bat
REM  (Belotti et al.'s bare PP-MISOCP is run separately from their
REM   GitHub repo; paste its numbers into the theirs_* columns.)
REM ============================================================

REM --- activate the conda environment that has gurobipy + the pipeline ---
REM     change "base" below if your environment has another name
call conda activate base

REM --- go to the code folder (edit if your path differs) ---
cd /d "C:\Users\hvalidi\Downloads\Sam_districting_paper\code"

REM --- run the benchmark (county + small tract instances, 1-hour limit each) ---
python run_A1_ours.py

echo.
echo Done. Results in A1_benchmark_results.xlsx ; maps in the maps\ folder.
pause
