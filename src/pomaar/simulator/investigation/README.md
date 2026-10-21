# Bistaticity investigation with idealized parametric beams

SBR+ study of how canonical targets' polarimetric signatures change between a monostatic and a
bistatic Tx/Rx geometry, i.e. what separating Tx and Rx phase centres in an automotive MIMO
array does to polarimetry (task note "Bistaticity analysis"; paper topic 4.1). 2026-10-20 to
2026-10-21; 600+ SBR+ solves. This folder holds the scripts and the documents:

| File | What |
| --- | --- |
| `investigation_log.md` | plan, dated log, distillation (findings, limitations, decisions), glossary (bistatic pattern, lobes, PTD/UTD, descriptors) |
| `rx_aim_reference.md` | the argument that receiver aiming changes only the span |
| `project.py` | the AEDT project and the data folders (below) |
| `build_scene.py`, `mesh.py`, `sweep.py`, `solve.py`, `export_sweep.py`, `reopen.py` | AEDT side: scenes, parametric sweeps, solves, exports |
| `run_production.sh`, `run_scaling.sh`, `run_extension.sh` | the runs behind the log |
| `summarize.py`, `pairs.py`, `scaling.py`, `aim_check.py`, `trust_checks.py` | host side (NumPy only): summaries, split pairs, scaling, receiver aiming, trust checks |
| `validation/` | the pipeline validation harness (separate project copy; its own README) |

## Scope: what "parametric beams" leaves out

Every antenna is an SBR+ parametric beam: 60° x 60°, identical H and V patterns, no pattern
cross-pol, a fixed phase centre, no mutual coupling and no array environment. These results
therefore describe the target and the geometry alone. They do **not** include:

* embedded element patterns, and the antenna diversity across a real array (each element's
  own pattern, gain and phase centre);
* the E-/H-plane gain difference between H and V elements, pattern cross-pol, and the rotation
  of the H/V basis off the principal planes, i.e. the receive polarization matrix P_Rx(β) of a
  real antenna;
* calibration errors, ground, multipath and clutter.

That isolation is also why the results are worth having. Whatever changes here is a property
of the target and the bistatic geometry, and real antennas add their own effects on top of it
rather than removing it. Statements of the form "X does not matter" (receiver aiming, splitting
the co-pol channels over separate pairs) hold only for these ideal beams. Each finding's
limitations are listed in the log's "Limitations" section.

## Main findings

* **Lobe-limited.** Strong-channel descriptors (Z_DR, ∠(S_HH S_VV*), entropy, ᾱ) hold up to
  0.6–1.0 β₀, where β₀ ≈ λ/D is the first null of the target's bistatic pattern. They change
  together at that null. The co-pol signature survives an array baseline d ≲ 0.5 Rλ/D. β₀
  scales as 1/D and is independent of range.
* **Weak channels move first.** The LDR of a target whose symmetry allows cross-pol (trihedral
  at roll 0° or 45°) rises at once: −39 dB at 1° for a 40 mm trihedral, falling about 9 dB per
  doubling of size. A real antenna's own cross-pol floor would mostly hide this.
* **Symmetry sets the order.** A mirror plane perpendicular to the displacement, or C2
  symmetry about the line of sight, makes the co-pol change only to second order in β. A
  mirror plane in the bistatic plane forbids cross-pol. A target with neither (slant
  trihedral) changes both co-pol and cross-pol to first order.
* **Receiver aiming** changes only the span (a scalar fit leaves at most −48 dB inside the
  lobe).
* **Z_DR from split pairs** (TxV–RxV, then TxH–RxH, the supervisor's layout) agrees with a
  co-located dual-polarized pair to within 0.2 dB and 0.7° for β ≤ 2°. It departs at larger
  angles, mostly for targets without the mirror or C2 symmetry above.

## Where everything lives

| What | Where | Tracked |
| --- | --- | --- |
| AEDT project (16 designs, 18 GB with results) | `~/Projects/AEDT/bistaticity-parametric-beams/bistaticity_parametric_beams.aedt` (`project.PROJECT_PATH`; `--project` takes another) | no |
| Exports (.npz) and run logs (`logs/`) | `data/interim/bistaticity/parametric_beams/` (`EXPORTS_DIR`) | no |
| Results: `summary/`, `pairs/`, `scaling/`, `rx_aim/` (figures and CSVs), `per_export/` (per-export analyses) | `data/processed/bistaticity/parametric_beams/` (`RESULTS_DIR`) | no |

The `parametric_beams` layer sits next to the full-wave cluster results in `data/*/bistaticity/`.
The results are regenerable from the exports with the host-side scripts.

| Results folder | What | Script |
| --- | --- | --- |
| `summary/` | signatures against β per target family (`dihedral`, `trihedral`, `sphere_plate`), receivers aimed (`aim1`) or on boresight (`aim0`); change from monostatic; aim difference; departure tables (CSV) | `summarize.py --aim 1`, `--aim 0` |
| `pairs/` | Z_DR and ∠(S_HH S_VV*) from split pairs against a co-located pair | `pairs.py` |
| `scaling/` | trihedral 40–160 mm at 1–16 m, dihedral 160 mm: β₀ and cross-pol against β and β/β₀ | `scaling.py <stems>` (stems in the log) |
| `rx_aim/` | the scalar fit behind `rx_aim_reference.md` | `aim_check.py` |
| `per_export/` | `bistaticity_analysis --per-receiver` of every `MultiRx*` export | see below |

Figures show unbiased quantities only (no reference lines or thresholds). Interpretive levels
appear only in the departure CSVs.

## Reproducing

Run everything from the repository root.

* **AEDT side** (in the container): `run_production.sh`, `run_scaling.sh` and
  `run_extension.sh` build, sweep, solve and export.
* **Host side**: `python src/pomaar/simulator/investigation/<script>.py` rewrites the results
  folders.
* **Per-export analyses**: run
  `bistaticity_analysis <export>.npz --per-receiver --out data/processed/bistaticity/parametric_beams/per_export`
  for each export.
