# Pipeline validation (2026-10-21)

End-to-end check of `hfss_array_builder` -> `sbr_simulator` -> `sbr_results` ->
`polarimetry_processor` / `bistaticity_analysis` against the running AEDT 2025 R2 GUI session,
on a copy of `mimo_polarimetry.aedt` (`data/interim/sandbox/validation/validation.aedt`, with
the results of `ApertureCoupledPatch` and the solved `ApertureCoupledPatchPolarimetricCluster`
only). `common.py` refuses any other project. The unit cell's sweep in the copy is cut to 5
points (`prepare.py`), so that a fresh array solves in minutes; the existing cluster keeps its
81-point solution.

A structured selection rather than every combination: each CLI option and each target type is
exercised at least once, each in a combination that changes something observable, plus the
error paths most likely to be hit.

## Plan and results

| ID | Stage | Case (what it exercises) | Result |
| --- | --- | --- | --- |
| B1 | builder | cluster layout, `--output-design`, `--use-existing-cs`, `--no-simulate` | Pass: placement verified, design validation passed, `Setup : Sweep` discrete 77–81 GHz (5 points) with radiated fields, port table written. |
| B2 | builder | same design again with `--overwrite --simulate` (replace path, solve, setup read-back) | Pass: the design is cleared in place and rebuilt; solved in 6 min 14 s. |
| B3 | builder | 3 x 4 contiguous ULA (7 V elements at yaw 0 and 180), `--no-simulate` | Pass: yaw-180 Tx ports get sign -1. The shipped `polarimetric_ula` layout is refused before any change (its Rx are `hv` and need a dual-port unit cell, `--extra-unit-cell hv=...`, which this project lacks). |
| B4 | builder | setup renamed to `ArraySetup`; a fresh builder process: `array_solution`, reports, `export_coupling_s_parameters` | Pass: `ArraySetup : Sweep` found and used for the reports and the Touchstone export. |
| B5 | builder | error: layout bound to a unit cell that the project lacks | **Bug, fixed** (921197b): PyAEDT's `Hfss(design=...)` created an empty design of that name, then the phase-centre fit crashed on an AEDT error. Now a `LayoutError` names the missing cell before anything is opened. |
| S1 | SBR+ | 5-point array, dihedral, x boresight (defaults), `--design` alias, `--no-simulate` | Pass: AntennaCS Euler (0, 90, 90), TargetCS expressions as in `scene_frames`, linked `ArraySetup : Sweep` (the renamed setup, found by itself), SBR+ sweep 77–81 GHz in 5 points, reports from the scattered solution. |
| S2 | SBR+ | z boresight, sphere, `--distance` alias, `--map-variables`, `--no-simulate` | Pass with `copolarSpacing` (declared as a `VariableProp` of the link). Mapping the dependent `clusterRadius` was refused correctly, but only after the design, CSs and target were created: **fixed** (e75f84d), now checked first, with an `[ERROR]` line instead of a traceback. |
| S3 | SBR+ | x boresight, dihedral roll 45 at (az 10, el 5), `--simulate`, export, steering | Pass: steering at (10, 5) removes 21.7 deg of geometric cross-pol phase (predicted from the virtual positions: 21.7 deg); residual +3.7 deg is antenna/target physics (patch patterns 10 deg off axis). |
| S4 | SBR+ | z boresight, same direction given as the scene's spherical angles (26.74, 11.17) | Pass with the processor's (10, 5): +3.6 deg, same as S3. Steering with the scene's own angles leaves +44.4 deg (predicted 41 deg off). **Fixed by design** (9359398): exports now record `boresight` and `target_steering_deg` (the processor's angles), `scene_frames.target_steering_angles` converts. |
| S5 | SBR+ | `--scene-only` on S1: trihedral at 3 m, roll 30 | Pass (with `--overwrite`; without it the prompt defaults to keeping the scene). A z design rebuilt without `--boresight z` would have got an x-style TargetCS: **fixed** (9359398), refused against AntennaCS. |
| S6 | SBR+ | 81-point solved cluster, no target (background), `--simulate`, `--overwrite` | |
| S7 | SBR+ | linking an unsolved array (ValUla) | Linked silently; the first SBR+ solve would have started the full-wave solve of the array ('simulate source design as needed'). Now warns (`warn_unsolved_source`). |
| S8 | SBR+ | error: a source design that does not exist | **Bug, fixed** (c61e9fe): it was created empty and then reported as having no linkable sweep. |
| E1 | export | `sbr_results` default (scattered) and `--setup 'Setup : Sweep'` (total), port table attached | Pass: port table found next to the project (signs corrected, positions known). Missing design or a project not open by name: clear `[ERROR]` (c61e9fe; before, PyAEDT created the missing design). |
| E2 | analysis | parametric sweep (targetRoll) -> `sbr_results --sweep-variable` -> `bistaticity_analysis` | |
