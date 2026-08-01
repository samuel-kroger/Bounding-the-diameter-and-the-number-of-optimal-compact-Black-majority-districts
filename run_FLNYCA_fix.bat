@echo off
setlocal EnableDelayedExpansion
REM ============================================================================
REM  Recovery run: redo ONLY what failed in the overnight sweep.
REM
REM      run_FLNYCA_fix.bat > flnyca_fix.log 2>&1
REM
REM  Finishes with  ===== FIX RUN COMPLETE =====   and does not pause.
REM
REM  WHAT ALREADY SUCCEEDED and is NOT repeated (about 10 hours of solver time):
REM    * instance info, both groups, except the CA majority-Latino row
REM    * ell_s for all six instances  (FL 11, NY 13, CA 11)
REM    * k_m for FL and NY, both groups, and CA majority-Black (k_m = 0)
REM    * the three-phase fixing tables
REM
REM  WHAT FAILED, WHY, AND WHAT WAS FIXED
REM  ------------------------------------
REM  1. CA majority-Latino k_m crashed with
REM         AttributeError: 'problem_instance' object has no attribute 'output'
REM     CA is the first instance ever to hit the one-hour limit on that model,
REM     and the TIME_LIMIT branch never set self.output. It also stored a raw
REM     float (33.0), which then broke range() on every later step. Both are
REM     fixed, and the timeout branch now reports the solver BOUND rather than
REM     the incumbent, which is the mathematically valid upper bound.
REM
REM  2. All six GerryChain runs failed to seed. The subprocess time-boxing did
REM     not work on Windows: every attempt ran to the 180 s timeout and was
REM     killed, burning 2 hours per instance for nothing, even though the same
REM     settings seed FL in about 5 s in-process. Seeding is now done in
REM     process by default. Verified: FL 5.7 s, NY 7.2 s, CA 15.2 s.
REM ============================================================================

call conda activate base
cd /d "%~dp0"

if "%GC_ITERS%"=="" set GC_ITERS=10000
set GC_SEED_SUBPROC=0
set GC_SEED_RETRIES=40
set GC_SEED_ATTEMPTS=20000

echo ############################################################
echo #  A  Sanity: seeding works before committing hours to it
echo ############################################################
python gc_seed.py FL 28
if errorlevel 1 goto :failed
python gc_seed.py NY 26
if errorlevel 1 goto :failed
python gc_seed.py CA 52
if errorlevel 1 goto :failed

echo.
echo ############################################################
echo #  B  Recover the CA majority-Latino k_m row  (~1.5 h)
echo ############################################################
python run_flny_table.py --table upper_bound_mm --group hispanic

echo.
echo ############################################################
echo #  C  Recover the CA row in the instance table  (fast, cached)
echo ############################################################
python run_flny_table.py --table instance_info --group hispanic

echo.
echo ############################################################
echo #  D  GerryChain, six instances in parallel, GC_ITERS=%GC_ITERS%
echo #     Finished instances are skipped, so this resumes cleanly.
echo ############################################################
for %%S in (FL NY CA) do (
  for %%G in (black hispanic) do (
    echo   launching %%S %%G
    start "gcfix_%%S_%%G" /min cmd /c "python run_gc_one.py --state %%S --group %%G > gc_log_%%S_%%G.txt 2>&1"
  )
)

echo.
echo   Waiting. Watch one with:  type gc_log_FL_hispanic.txt
:waitloop
timeout /t 60 /nobreak > nul
tasklist /fi "windowtitle eq gcfix_*" 2>nul | find /i "cmd.exe" > nul
if not errorlevel 1 goto :waitloop

python merge_gc.py

echo.
echo ############################################################
echo #  E  Frontier, now that incumbents exist to warm start it
echo ############################################################
python run_flny_table.py --table point_check --group black
python run_flny_table.py --table point_check --group hispanic

echo.
echo ############################################################
echo #  RESULTS
echo ############################################################
type gerrychain_FLNYCA.csv 2>nul
echo.
dir /b results_black\*_FLNYCA_*.csv    2>nul
dir /b results_hispanic\*_FLNYCA_*.csv 2>nul
echo.
echo ===== FIX RUN COMPLETE =====
exit /b 0

:failed
echo.
echo *** Seeding sanity check FAILED. Stopping before any solver work. ***
echo ===== FIX RUN FAILED =====
exit /b 1
