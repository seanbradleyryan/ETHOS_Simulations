"""
calc_beam_plan_doses_allsessions.py  -  run calc_beam_plan_doses.py on every session case

PURPOSE
    Each fraction is imported as its own case ("Session_1" ... "Session_10",
    see import_raystation_input.py), so calc_beam_plan_doses.py had to be
    started by hand once per case. This script makes each "Session_N" case
    current in turn (case.SetCurrent) and runs calc_beam_plan_doses.py on it
    unchanged, so the per-case calculations are identical.

EXECUTION (Raystation 2024 only).
    Run inside RayStation (Scripting > Script creation > Run) with the
    patient open; any case may be current. It can be resumed:
    calc_beam_plan_doses.py skips plan/CT pairs already DONE in each
    session's progress log, so a re-run only does the unfinished work.

INPUTS
    - Cases named "Session_N" (other cases are ignored). No manual case
      preparation: calc_beam_plan_doses.py runs prepare_case (Generic CT,
      duplicate ROIs, material overrides, external ROI) on each case before
      its dose calculation.
    - CALC_SCRIPT_PATH below: location of calc_beam_plan_doses.py.

OUTPUTS
    Same as calc_beam_plan_doses.py, for every session case.
    calc_beam_plan_doses.py prints its own errors and returns, so an error
    in one session does not stop the remaining sessions. A case that cannot
    be made current is skipped and listed at the end.
"""

from connect import *
import os
import re
import runpy


# ============================================================
# Configuration
# ============================================================
CALC_SCRIPT_PATH = "F:/ETHOS_Simulations/RayStation/calc_beam_plan_doses.py"

CASE_NAME_RE = re.compile(r'^Session_(\d+)$', re.IGNORECASE)


# ============================================================
# Main
# ============================================================
if not os.path.isfile(CALC_SCRIPT_PATH):
    raise RuntimeError(f"calc_beam_plan_doses.py not found at {CALC_SCRIPT_PATH}")

patient = get_current("Patient")

# Session cases, sorted by number (Session_2 before Session_10)
session_cases = []
for case in patient.Cases:
    m = CASE_NAME_RE.match(case.CaseName)
    if m:
        session_cases.append((int(m.group(1)), case))
session_cases.sort(key=lambda x: x[0])

if not session_cases:
    raise RuntimeError("No cases named 'Session_N' found for this patient.")

print(f"Found {len(session_cases)} session case(s): "
      + ", ".join(c.CaseName for _, c in session_cases))

skipped_cases = []
for _, case in session_cases:
    case_name = case.CaseName
    print(f"\n{'#'*60}")
    print(f"CASE: {case_name}")
    print(f"{'#'*60}")

    try:
        case.SetCurrent()
        current_name = get_current("Case").CaseName
        if current_name != case_name:
            raise RuntimeError(f"current case is still '{current_name}'")
    except Exception as exc:
        print(f"WARNING: could not open case '{case_name}': {exc}. Skipping.")
        skipped_cases.append(case_name)
        continue

    # calc_beam_plan_doses.py reads get_current("Case"), i.e. the case set above
    runpy.run_path(CALC_SCRIPT_PATH, run_name="__main__")

print(f"\n{'#'*60}")
print(f"All session cases processed. Skipped: {skipped_cases if skipped_cases else 'none'}")
print(f"{'#'*60}")
