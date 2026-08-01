@echo off
REM ===================================================================
REM  Table 6 / Figure -- the (s, #MM) tradeoff frontier
REM
REM  For each instance and each k_m, finds s*(k_m) = the smallest diameter
REM  bound at which a plan with k_m majority-minority districts exists.
REM  (Reviewer 1 Comment 3 / Reviewer 2 Major 3.)
REM
REM  Uses all three accelerations together:
REM    * heuristic (GerryChain) -- a VERIFIED plan certifies feasibility with
REM      no MIP at all (population balance, diameter <= s, and >= k_m MM
REM      districts are each checked explicitly before the plan is accepted);
REM    * variable fixing        -- shrinks the prescribed feasibility MIP;
REM    * symmetry breaking      -- what makes the INFEASIBLE s values provable
REM      (the regime where symmetry pays off, cf. Table 5).
REM
REM  Monotonicity: a plan with k_m MM districts is also a plan for k_m - 1, so
REM  s*(k_m) is non-decreasing in k_m. The sweep runs k_m upward and carries s
REM  forward, skipping every s already known infeasible. It is capped at u_s,
REM  so it always terminates.
REM
REM  STEP 1 builds the heuristic plans (pickles under results_<group>\gerrychain\).
REM  STEP 2 computes the frontier.  Output: results_<group>\point_check.csv
REM  Resumable: rows are keyed on (state, level, k_m); re-running continues.
REM ===================================================================
call conda activate base
cd /d "%~dp0"

echo === STEP 1: GerryChain heuristic plans (feasibility certificates) ===
python run_table.py --table gerrychain --group black    --dataset paper
python run_table.py --table gerrychain --group hispanic --dataset hispanic_paper

echo.
echo === STEP 2: (s, #MM) frontier ===
python run_table.py --table point_check --group black    --dataset paper
python run_table.py --table point_check --group hispanic --dataset hispanic_paper

echo.
echo Wrote results_black\point_check.csv and results_hispanic\point_check.csv
echo (Plot the (s, #MM) frontier from these CSVs for the tradeoff figure.)
pause
