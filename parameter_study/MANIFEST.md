# Parameter Study Manifest

Background work completed before the Task 1/2/3 assignment from Yongsheng.
Comparing solitary vs regular waves on a monopile/cylinder across wave heights,
periods, and wave theories. Presented to manager prior to current assignment.

| Case | Week | Status | Notes |
|---|---|---|---|
| solitary3D | Week 1 (Jun 16) | Superseded | Early exploratory case, replaced by solitary_H* sweep |
| solitary3DPorous | Week 1 (Jun 15) | Superseded | Early exploratory case (porous variant), replaced by solitary_H* sweep |
| monopileSolitary | Week 2 (Jun 24) | Completed | |
| stokes2_T3 | Week 2 (Jun 22) | Completed | |
| regular_T2 | Week 3 (Jun 22) | Completed | Period sweep |
| regular_T3 | Week 3 (Jun 22) | Completed | Period sweep |
| regular_T5 | Week 3 (Jun 22) | Completed | Period sweep |
| regular_T8 | Week 3 (Jun 22) | Completed | Period sweep |
| solitary_H005 | Week 3 (Jun 22) | Completed | Wave height sweep |
| solitary_H015 | Week 3 (Jun 22) | Completed | Wave height sweep |
| solitary_H020 | Week 3 (Jun 22) | Completed | Wave height sweep |
| regular_T3_H005 | Week 3 (Jun 24) | Completed | Wave height sweep at T3 |
| regular_T3_H015 | Week 3 (Jun 24) | Completed | Wave height sweep at T3 |
| regular_T3_H020 | Week 3 (Jun 24) | Completed | Wave height sweep at T3 |
| chappelear_H005 | Week 3 (Jun 26) | Scrapped | Anomalously high forces (1100-1150N), briefly investigated, abandoned as not priority |
| chappelear_H010 | Week 3 (Jun 26) | Scrapped | Same as above |
| chappelear_H015 | Week 3 (Jun 26) | Scrapped | Same as above |
| chappelear_H020 | Week 3 (Jun 26) | Scrapped | Same as above |

## Key finding
Solitary vs regular wave force crossover at H≈0.15m.

## Scripts
- `compare_cases.py` — compares all cases by reading postProcessing/forces
- `extract_eta.py`, `plot_results.py` — elevation/result plotting utilities
- `run_sweep.sh` — automation for running the parameter sweep
- `sweep_results.txt` — raw sweep output
- `figures/` — result plots presented to manager
