@echo off
REM ============================================================================
REM  Fill the enacted-plan columns for FL, NY and CA.
REM
REM      run_enacted_FLNYCA.bat > enacted.log 2>&1
REM
REM  Ends with  ===== ENACTED COMPLETE =====  and does not pause.
REM  Needs internet: it downloads the 118th Congress district shapefiles from
REM  the Census TIGER site. Everything else is local and cached.
REM
REM  WHY: Table 3's last three columns (enacted diameter s, split parcels,
REM  enacted #MM) come from raw_data/current_plan/<STATE>.shp, which exists for
REM  the 28 states already in the paper but not for FL, NY or CA. Those three
REM  currently read NA, and since every other row in the table is populated,
REM  the gap would stand out. Total run time is a few minutes; nothing is
REM  re-solved, because ell_s and k_m are read from cache.
REM ============================================================================

call conda activate base
cd /d "%~dp0"

echo ############################################################
echo #  1  Download the enacted 118th Congress plans
echo ############################################################
python fetch_current_plan.py FL NY CA
if errorlevel 1 goto :failed

echo.
echo   Shapefiles now present:
dir /b raw_data\current_plan\FL.shp raw_data\current_plan\NY.shp raw_data\current_plan\CA.shp 2>nul

echo.
echo ############################################################
echo #  2  Rebuild the instance table so the columns populate
echo ############################################################
python run_flny_table.py --table instance_info --group black
python run_flny_table.py --table instance_info --group hispanic

echo.
echo ############################################################
echo #  3  Result
echo ############################################################
type results_black\instance_info_FLNYCA_black.csv
echo.
type results_hispanic\instance_info_FLNYCA_hispanic.csv
echo.
echo If the last three columns still read NA, the download did not land;
echo check enacted.log for the URL that failed.
echo.
echo ===== ENACTED COMPLETE =====
exit /b 0

:failed
echo.
echo *** Download failed. Check your connection and the log. ***
echo ===== ENACTED FAILED =====
exit /b 1
