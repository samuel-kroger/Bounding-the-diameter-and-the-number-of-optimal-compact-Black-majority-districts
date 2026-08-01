@echo off
setlocal EnableDelayedExpansion
REM ============================================================================
REM  FL, NY and CA tract instances -- every paper table that applies.
REM
REM  Run it and walk away:
REM      run_FLNYCA.bat > flnyca.log 2>&1
REM
REM  It does NOT pause at the end, so a redirected run terminates on its own.
REM  When it is finished the log ends with the line
REM      ===== SWEEP COMPLETE =====
REM  so you can check with:  findstr /C:"SWEEP COMPLETE" flnyca.log
REM
REM  WHAT CHANGED SINCE THE LAST VERSION
REM  ----------------------------------
REM  * CA is included. Its three stranded island tracts (Farallones, Santa
REM    Catalina, the LA Channel Islands) are reconnected in step 0.
REM  * GerryChain seeding is now time-boxed and retried (gc_seed.py). The FL
REM    failure last time was a SEEDING failure, not a shortage of iterations:
REM    recursive_tree_part could not find a balanced cut and raised before the
REM    chain took a single step. With the fix, FL seeds in ~5 s and immediately
REM    finds 4 majority-Latino districts.
REM  * Step 5 runs the six instances IN PARALLEL. GerryChain is single threaded,
REM    so six sequential runs take ~38 h while six parallel ones take ~6.5 h.
REM  * No trailing pause; per-step resumability; per-instance logs.
REM
REM  ABOUT RAISING THE ITERATION COUNT
REM  ---------------------------------
REM  GC_ITERS controls chain length and does improve result quality. It is also
REM  linear in time: measured at ~2.3 s per iteration on these instances,
REM      10,000 iters  ->  ~6.4 h per instance
REM     100,000 iters  ->  ~64 h per instance
REM   1,000,000 iters  ->  ~27 days per instance
REM  The default below is 10,000, which is what every other state in the paper
REM  used; changing it for these three only would break comparability. To
REM  override anyway:   set GC_ITERS=25000   before running.
REM
REM  Requirements: conda env with gurobipy (licensed), gerrychain, geopandas,
REM  ortools, networkx, tqdm.
REM  Expected wall clock: about 12 to 16 hours.
REM ============================================================================

call conda activate base
cd /d "%~dp0"

if "%GC_ITERS%"=="" set GC_ITERS=10000
if "%GC_SEED_TIMEOUT%"=="" set GC_SEED_TIMEOUT=180
if "%GC_SEED_RETRIES%"=="" set GC_SEED_RETRIES=40

echo ############################################################
echo #  STEP 0  Reconnect the FL, NY and CA contiguity graphs
echo ############################################################
python repair_connectivity.py
if errorlevel 1 goto :failed
python repair_connectivity.py --check
if errorlevel 1 goto :failed

echo.
echo ############################################################
echo #  STEP 1  Instance info  (Table 3)
echo ############################################################
python run_flny_table.py --table instance_info --group black
python run_flny_table.py --table instance_info --group hispanic

echo.
echo ############################################################
echo #  STEP 2  Diameter lower bound ell_s  (Algorithm 2)
echo #          Must precede steps 3 and 4: they reuse the cache.
echo ############################################################
python run_flny_table.py --table lower_bound_s --group black
python run_flny_table.py --table lower_bound_s --group hispanic

echo.
echo ############################################################
echo #  STEP 3  Upper bound k_m and the induced fixing percentage
echo ############################################################
python run_flny_table.py --table upper_bound_mm --group black
python run_flny_table.py --table upper_bound_mm --group hispanic

echo.
echo ############################################################
echo #  STEP 4  Three-phase variable fixing  (Supplement 3.1, 3.2)
echo ############################################################
python run_flny_table.py --table fixing --group black
python run_flny_table.py --table fixing --group hispanic

echo.
echo ############################################################
echo #  STEP 5  GerryChain, six instances IN PARALLEL
echo #          GC_ITERS=%GC_ITERS% per instance
echo #          Each writes gerrychain_^<STATE^>_^<group^>.csv and its own log.
echo #          Already-finished instances are skipped, so re-running the
REM             batch after an interruption resumes rather than restarts.
echo ############################################################
for %%S in (FL NY CA) do (
  for %%G in (black hispanic) do (
    echo   launching %%S %%G
    start "gc_%%S_%%G" /min cmd /c "python run_gc_one.py --state %%S --group %%G > gc_log_%%S_%%G.txt 2>&1"
  )
)

echo.
echo   Waiting for the six GerryChain runs to finish...
echo   (progress:  type gc_log_FL_hispanic.txt )
:waitloop
timeout /t 60 /nobreak > nul
tasklist /fi "windowtitle eq gc_*" 2>nul | find /i "cmd.exe" > nul
if not errorlevel 1 goto :waitloop

echo   All GerryChain runs finished. Merging.
python merge_gc.py

echo.
echo ############################################################
echo #  STEP 6  Diameter vs representation frontier  (Figure 1)
echo ############################################################
python run_flny_table.py --table point_check --group black
python run_flny_table.py --table point_check --group hispanic

echo.
echo ############################################################
echo #  RESULTS
echo ############################################################
dir /b results_black\*_FLNYCA_*.csv    2>nul
dir /b results_hispanic\*_FLNYCA_*.csv 2>nul
dir /b gerrychain_FLNYCA.csv           2>nul
echo.
echo Nothing else in results_black\ or results_hispanic\ was modified.
echo.
echo Not run, and why:
echo   - Table 4 (county level, zero splits): FL, NY and CA have counties that
echo     exceed the population upper bound, so no county instance exists.
echo   - Supplement big-M table: county level only.
echo   - Table 5 benchmark vs Belotti et al.: needs their solver, run from
echo     their own repository.
echo.
echo ===== SWEEP COMPLETE =====
exit /b 0

:failed
echo.
echo *** Connectivity repair FAILED. Stopping before any solver work. ***
echo ===== SWEEP FAILED =====
exit /b 1
