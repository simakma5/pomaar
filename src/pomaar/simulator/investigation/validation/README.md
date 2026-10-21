# Pipeline validation (2026-10-21)

End-to-end check of `hfss_array_builder` -> `sbr_simulator` -> `sbr_results` ->
`polarimetry_processor` / `bistaticity_analysis` against the running AEDT 2025 R2 GUI session,
on a copy of `mimo_polarimetry.aedt`
(`~/Projects/AEDT/mimo-polarimetry-validation/mimo_polarimetry_validation.aedt`, with
the results of `ApertureCoupledPatch` and the solved `ApertureCoupledPatchPolarimetricCluster`
only). `common.py` refuses any other project. The unit cell's sweep in the copy is cut to 5
points (`prepare.py`), so that a fresh array solves in minutes; the existing cluster keeps its
81-point solution. The cases' logs and exports are in `data/interim/validation/` (untracked). The
X11 helpers (`x11_*.py`, `watch_dialogs.sh`, for AEDT's modal dialogs on the VNC display) need
the project's x11 extra: `uv sync --extra x11`.

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
| S6 | SBR+ | 81-point solved cluster, no target (background), `--simulate`, `--overwrite` (replaces the existing `LinkedApertureCoupledPatchPolarimetricCluster`) | **Two bugs, fixed.** SBR+ refused to solve ('There should be at least one model object in SBR+ design'), so `--target none` had never been solvable; it now gets a 1 mm vacuum cube at TargetCS (2b6b0c7): scattered solution -225 dB (nothing), total -20.5 dB co-pol / -25.4 dB cross-pol (the direct coupling), 81 frequencies. The failed solve exited with 0: failed solves now raise `SolveFailedError`, exit code 3, in the builder and `sbr_simulator` (df1e9f1). |
| S7 | SBR+ | linking an unsolved array (ValUla) | Linked silently; the first SBR+ solve would have started the full-wave solve of the array ('simulate source design as needed'). Now warns (`warn_unsolved_source`). |
| S9 | SBR+ | `--overwrite` of the same SBR+ design twice in a row | **Bug, fixed** (7cbf828): the second run hung on AEDT's modal 'The file <design>.asol_priv is in use - do you want to wait?' (renaming before deleting is not enough once the design has result files). `sbr_simulator` now saves, closes and reopens the project after deleting; two overwrites in a row pass. |
| S8 | SBR+ | error: a source design that does not exist | **Bug, fixed** (c61e9fe): it was created empty and then reported as having no linkable sweep. |
| E1 | export | `sbr_results` default (scattered) and `--setup 'Setup : Sweep'` (total), port table attached | Pass: port table found next to the project (signs corrected, positions known). Missing design or a project not open by name: clear `[ERROR]` (c61e9fe; before, PyAEDT created the missing design). |
| E2 | analysis | parametric sweep (targetRoll 0 ... 90 by 22.5) of S3 -> `sbr_results --sweep-variable` -> `bistaticity_analysis --azimuth 10 --elevation 5` | Pass: 5 x 4 x 5 export, figures and CSV; roll 0/90 co-pol with arg(S_HH S_VV*) ~ 180 deg and cross-pol 21 dB down, roll 45 cross-pol only, 22.5/67.5 all four channels equal. |

## Findings

Bugs (all fixed and re-tested):

1. PyAEDT's `Hfss(design=name)` creates an empty design when `name` does not exist. The builder
   (missing unit cell), `sbr_simulator` (missing source or `--scene-only` design) and
   `sbr_results` (missing design) each did so, then failed far from the cause. All now check
   with `require_designs` / `project_design_names` first.
2. `sbr_simulator --target none` could never be solved (no model object). Vacuum placeholder.
3. Failed solves exited with 0. Now `SolveFailedError` and exit code 3.
4. `--overwrite` of an SBR+ design that has result files hung on a modal dialog. The project is
   now reopened after deleting.
5. Invalid `--map-variables` were caught only after the design was built. Checked first now.
6. `--scene-only` on a z-boresight design without `--boresight z` built an x-style TargetCS.
   Refused against AntennaCS now.

Conventions clarified:

7. With a z boresight, the scene's targetAzimuth/targetElevation are spherical angles (azimuth
   about boresight, polar angle off it), not the processor's steering angles (towards H, towards
   V). S4 showed a 41 deg cross-pol phase error from mixing them. Exports now carry
   `target_steering_deg` in the processor's angles (and `boresight`).

Usability:

8. The layout generators' `-o` resolved relative paths against the generator's folder (1b3c211).
9. Linking an unsolved array now warns that the first SBR+ solve runs its full-wave solve.
10. The shipped `polarimetric_ula` layout needs a dual-port unit cell (`--extra-unit-cell
    hv=<Design>`), which `mimo_polarimetry.aedt` lacks; the builder says so before any change.

Not covered: `--field-type nearfield`, `--roughness`/`--height-deviation`, `phase_reference_mode:
single` layouts, dual-port (`hv`) unit cells, an interrupted solve.
