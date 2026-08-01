@echo off
REM ============================================================================
REM  FL and NY tract instances -- full sweep of every paper table that applies.
REM
REM  Context: CA, FL, NY were excluded from the tract-level experiments because
REM  their RAW contiguity graphs are disconnected. Step 0 repairs FL and NY by
REM  adding a handful of explicit water-crossing edges, following Validi et al.
REM  (2022). CA is deliberately left out: at n = 9,129 and k = 52 it is much
REM  larger than TX (n = 6,896, k = 38), our largest instance.
REM
REM  SAFETY: every table is run through run_flny_table.py, which snapshots the
REM  existing results CSVs, runs the FL/NY sweep, writes the new rows to
REM  <name>_FL_NY.csv, and restores the originals. The numbers behind the tables
REM  already in the manuscript are never overwritten.
REM
REM  Requirements: conda env with gurobipy (licensed), gerrychain, geopandas,
REM  ortools, networkx, tqdm.
REM
REM  Expected wall clock: roughly 12 to 30 hours. The one-hour solver limit
REM  applies per MIP, and several steps solve a sequence of MIPs. Run it
REM  overnight. Every step is safe to re-run; the caches are reused.
REM ============================================================================

call conda activate base
cd /d "%~dp0"

echo.
echo ############################################################
echo #  STEP 0  Repair the FL and NY contiguity graphs
echo ############################################################
python repair_connectivity.py
if errorlevel 1 goto :failed
python repair_connectivity.py --check
if errorlevel 1 goto :failed

echo.
echo ############################################################
echo #  STEP 1  Table 3, instance info (n, m, diameter, %% minority VAP)
echo ############################################################
python run_flny_table.py --table instance_info --group black
python run_flny_table.py --table instance_info --group hispanic

echo.
echo ############################################################
echo #  STEP 2  Diameter lower bound ell_s  (Algorithm 2)
echo #          Table 6 column 1. MUST run before steps 3 to 5:
echo #          they reuse the cached ell_s and independent sets.
echo ############################################################
python run_flny_table.py --table lower_bound_s --group black
python run_flny_table.py --table lower_bound_s --group hispanic

echo.
echo ############################################################
echo #  STEP 3  Upper bound k_m on majority-minority districts
echo #          plus the induced variable-fixing percentage.
echo #          Table 6 columns k_m and fix%%.
echo ############################################################
python run_flny_table.py --table upper_bound_mm --group black
python run_flny_table.py --table upper_bound_mm --group hispanic

echo.
echo ############################################################
echo #  STEP 4  Three-phase variable fixing (Online Supplement 3.1, 3.2)
echo #          Reports the strengthened bound ell_s(1) and %% fixed.
echo ############################################################
python run_flny_table.py --table fixing --group black
python run_flny_table.py --table fixing --group hispanic

echo.
echo ############################################################
echo #  STEP 5  GerryChain short bursts: incumbents and warm starts.
echo #          Feeds Table 6 (inc., s_inc) and Figure 1.
echo ############################################################
python run_flny_table.py --table gerrychain --group black
python run_flny_table.py --table gerrychain --group hispanic

echo.
echo ############################################################
echo #  STEP 6  Diameter vs representation frontier (Figure 1)
echo ############################################################
python run_flny_table.py --table point_check --group black
python run_flny_table.py --table point_check --group hispanic

echo.
echo ############################################################
echo #  DONE. New files, all suffixed _FL_NY:
echo ############################################################
dir /b results_black\*_FL_NY.csv    2>nul
dir /b results_hispanic\*_FL_NY.csv 2>nul
dir /b gerrychain_FL_NY.csv         2>nul
echo.
echo Nothing else in results_black\ or results_hispanic\ was modified.
echo.
echo NOT run here, and why:
echo   - Table 4 (county-level, zero splits): FL and NY have no feasible
echo     county-level instance at +/-0.5%% deviation.
echo   - Supplement big-M table: county-level only.
echo   - Supplement symmetry table: very expensive. To add it, run
echo       python run_flny_table.py --table symmetry --group hispanic
echo   - Table 5 benchmark vs Belotti et al.: needs their solver run
echo     separately from their repository.
goto :done

:failed
echo.
echo *** Connectivity repair FAILED. Stopping before any solver work. ***
pause
exit /b 1

:done
echo.
pause
