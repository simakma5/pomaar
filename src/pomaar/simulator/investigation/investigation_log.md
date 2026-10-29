# Bistaticity investigation: running plan and log

Working notes of the sandboxed investigation behind the task note "Bistaticity analysis"
(Obsidian). A distillation goes to that note's Findings and Log; this file keeps the details.
All scenes use idealized parametric beams (see `README.md` next to this file for the scope).
Results folders such as `summary/`, `pairs/`, `scaling/` and `rx_aim/` below are in
`data/processed/bistaticity/parametric_beams/`; until 2026-10-21 they lived in the sandbox data
folder (`<sandbox>/...` in older entries), which no longer exists.

## Sandbox

| What | Where |
| --- | --- |
| Branch | `sandbox/bistaticity`, worked in a git worktree `~/Repositories/pomaar-sandbox` (`main` and its working tree untouched); merged into `main` and the worktree removed on 2026-10-21. Since then on `main`, one topical commit per step (user: "work in the main branch, committing everything in a structured manner") |
| Baseline | the uncommitted simulator state of `main` on 2026-10-20, first one commit (`9dfe44b`), then split into six topical commits on `main` (`d1af11f` ... `7ac4186`) onto which the branch was rebased |
| AEDT project | `~/Projects/AEDT/bistaticity-parametric-beams/bistaticity_parametric_beams.aedt` (until 2026-10-21 `data/interim/sandbox/bistaticity/bistaticity_sandbox.aedt`), a copy of `~/Projects/AEDT/bistaticity-analysis/bistaticity_analysis.aedt` (12:21), open in the GUI session next to the original |
| Scripts | `src/pomaar/simulator/investigation/` (run with `ansys-vnc exec python ...` from the repository root); they refuse any project but the sandbox copy |
| Exports | `data/interim/bistaticity/parametric_beams/` (.npz, untracked; run logs in `logs/`) |
| Results | `data/processed/bistaticity/parametric_beams/` (untracked, regenerable): `summary/`, `pairs/`, `scaling/`, `rx_aim/`, and `per_export/<export>/` |
| Scripts and documents | `src/pomaar/simulator/investigation/` (tracked): this log, `README.md`, `rx_aim_reference.md` |

## Current action point

> Put multiple receivers in the scene (at least 3: monostatic, reasonably bistatic, and twice
> that).

### Scene design (parametric beams, x boresight, ISO 8855)

* One dual-polarized Tx site at the origin and three dual-polarized Rx sites, all ideal
  parametric beams (no coupling, no pattern cross-pol, no phase-centre wander), so that only
  the bistatic geometry and the target remain.
* Rx site k (k = 0, 1, 2) sits at y_k = targetRange tan(k bistaticAngle): for a boresight
  target, its bistatic angle is exactly k bistaticAngle (0 = monostatic, co-located with Tx).
  One solve gives the three receivers; sweeping `bistaticAngle` scans "reasonably bistatic".
  Receivers with equal k bistaticAngle from different solves (e.g. Rx_2 at 1 deg and Rx_1 at
  2 deg) sit at the same position (y depends on k bistaticAngle only), so their agreement only
  checks the bookkeeping of the export and the pooling, not the physics.
* `rxAim` (0 or 1) turns each receiver towards the target (yaw -k bistaticAngle), so that the
  antenna pattern and its off-axis polarization definition drop out (1), or keeps every
  antenna on boresight like a real array (0). The difference separates the geometric
  polarization squint (a coordinate artefact) from the bistatic scattering itself.
* `spacingHV` (H at z = +s/2, V at z = -s/2 in every site) stays 0 here (dual-pol): the
  polarization channels share each site's phase centre.
* Names: instances `Tx_H`, `Tx_V`, `Rx_<k>_H`, `Rx_<k>_V`; SBR+ appends `_p1` to their ports,
  which `parse_port` now strips.

### Plan

* [x] M0 Sandbox: worktree and branch, project copy, plan and log.
* [x] M1 Code: `parse_port` accepts SBR+ native-antenna ports (`<component><index>_p<N>`);
  sandbox helpers (`investigation/`).
* [x] M2 Build the multi-receiver scene in the sandbox project (`build_scene.py`;
  `MultiRxTrihedral`: 40 mm trihedral at 4 m, bistaticAngle 1 deg, rxAim 1, spacingHV 0).
* [x] M3 Trust checks (small solves), then **stop for validation** (all pass; see the log):
  * T1 level: monostatic trihedral |S_HH| = |S_VV| against the radar equation;
  * T2 polarization labels: relative to a trihedral, a dihedral with its seam along H flips
    S_VV and keeps S_HH (the planned PEC strip along H turned out uninformative);
  * T3 phase: arg S follows -k (r_Tx + r_Rx) of each antenna's own position;
  * T4 numerical cross-pol floor against ray density (4, 8, 16 rays per wavelength);
  * T5 the co-located monostatic receiver (Rx_0 on top of Tx) solves cleanly;
  * T4b symmetry: the trihedral's cross-pol in the bistatic plane vanishes when that plane is
    its mirror plane (targetRoll 90 deg) and may not otherwise (targetRoll 0).
* [x] M4 Production sweep (validated 2026-10-20, with a dihedral at 22.5 deg added):
  * designs (all 40 mm across, at 4 m): `MultiRxSphere` (radius 20 mm, mesh slider 9, T6),
    `MultiRxTrihedral` at targetRoll 0 and 90 deg, `MultiRxDihedral` at targetRoll 0, 22.5,
    45, 67.5 and 90 deg (a linear roll sweep; 67.5 mirrors 22.5);
  * `bistaticAngle` 0.25 ... 8 deg in 0.25 deg steps (32 values; with Rx_2 the bistatic
    angles reach 16 deg), each with rxAim 0 and 1: 512 solves, about 3.5 h
    (`investigation/run_production.sh`, log `logs/production.log` next to the exports);
  * exports: one per (design, rxAim, targetRoll) over `bistaticAngle`, with the scene
    geometry (`positions_m`, `target_m`) in the metadata (`export_sweep.py`);
  * analysis: `bistaticity_analysis --per-receiver` pools every (receiver, variation) by its
    bistatic angle (`pool_receivers`); a summary script compares targets and rxAim.
  * [x] T6 sphere faceting: at the default mesh (slider 5) the 20 mm sphere was 2.4 dB above
    the radar equation with +-2.5 dB ripple over the band, and its co-pol ratio drifted
    0.18 dB / 1.7 deg at 4 deg bistatic; at slider 8 and 9: -113.1 dB against -113.3 dB,
    +-0.6 dB ripple, 0.027 dB / 0.43 deg at 4 deg (8 and 9 agree), cross-pol floor -90 / -95
    dB. The sphere runs at slider 9; the plate targets are not faceted by curvature.
  * [x] T7 symmetry in the sign of the bistatic angle (receivers mirrored to -y, 3 solves): a
    target with C2 symmetry about the line of sight gives S(-beta) = S(beta), one with a
    mirror plane x-z gives S(-beta) = M S(beta) M with M = diag(-1, 1) (H flips).
* [x] M5 Summary of the action point for the task note (Findings, Log): see "Distillation".
* [x] M6 Size and range scaling (user request after M5): trihedral 40 / 80 / 160 mm at 1 / 4 /
  16 m, dihedral 160 mm; T8 edge-diffraction model (PO, PTD, PTD + UTD).
* [x] M7 Pipeline validation (user batch of 2026-10-21, task 1): builder -> SBR+ -> export ->
  analysis on a copy of `mimo_polarimetry.aedt`, 16 cases, 6 bugs fixed; plan and results in
  `investigation/validation/README.md`.
* [x] M8 Supervisor's suggestions (task 2): LDR and cross-co phase in radar conventions, more
  reflectors (slant trihedral at roll 45, flat plate), and Z_DR from split antenna pairs
  (`build_scene --layout pairs`, `pairs.py`); figures split by signature group.
* [x] M9 Consolidation (task 3): decisions, take-aways and limitations, glossary (below).

### Extension action points (not in this round)

* [ ] Full-wave check of the trihedral's bistatic cross-pol (HFSS IE / MoM, 40 mm edges at
  79 GHz, about 10 lambda): PO, PTD and UTD differ by 10 dB.
* [ ] Vertical receiver displacement (bistatic plane x-z, along V), next to the horizontal
  one (x-y, along H) studied here: the trihedral's and dihedral's symmetry planes swap roles.
* [ ] Real antenna patterns instead of parametric beams (the array builder's linked far fields):
  E/H-plane gain difference, pattern cross-pol and basis rotation off the principal planes,
  i.e. the receive polarization matrix P_Rx(beta) that the ideal beams leave out.
* [ ] Split pairs with real, non-identical H and V antennas and a per-pair calibration, the case
  where the split might matter in practice (the ideal beams cancel by mirror symmetry).
* [ ] Further reflectors if needed: a cylinder (singly curved, like a pole) and a thin wire or
  dipole (strongly polarizing, cross-pol at 45 deg); not run.

## Log

### 2026-10-20

* M0: created the worktree and branch `sandbox/bistaticity` from `main` (2afa6be) and
  committed the uncommitted simulator state of `main` as the baseline (9dfe44b). Copied the
  project (saved 12:21) to the sandbox as `bistaticity_sandbox.aedt`: AEDT cannot open two
  projects of the same name, and the distinct name keeps every script away from the original.
* M1/M2: AEDT names a parametric-beam instance '<component>1' and its port '<instance>_p1';
  PyAEDT lists excitations with a mode suffix ':1', which `SetSBRTxRxSettings` silently
  ignores (empty Tx/Rx matrix). `parse_port` and the builder handle both.
* Pitfall: re-creating a design under a just-deleted name (even after renaming it away) hits
  a modal "file <design>.asol_priv is in use" dialog on save, because AEDT keeps the deleted
  design's result files open. Dismissed it with XTEST clicks (python-xlib in the sandbox venv,
  helper scripts in the sandbox data folder); the builder now saves, closes and reopens the
  sandbox project before reusing a name (`reopen_sandbox`).
* A nominal solve of the 3-receiver scene (2 Tx x 6 Rx ports, 1001 frequencies, 4 rays per
  wavelength) takes 25 s: large sweeps are cheap.
* T1 level (`t1_trihedral_1deg`, Rx_0): |S_HH| = |S_VV| = -86.0 ... -84.8 dB over the band,
  mean -85.3 dB, against -85.6 dB from the radar equation (sigma = 4 pi a^4 / (3 lambda^2) =
  -1.3 dBsm at 79 GHz, G = 41253 / (60 x 60) = 10.6 dBi): +0.25 dB. Pass.
* T3 paths (phase slope over 74-84 GHz): 7999.94 / 8000.56 / 8002.42 mm for Rx_0/1/2 against
  r_t + r_r to the trihedral corner of 8000.00 / 8000.61 / 8002.44 mm (within 0.06 mm), the
  same for HH and VV. Every antenna's phase follows its own position: SBR+'s
  'SinglePhaseCenter' mode does not collapse them. Pass. (The image-through-the-corner path
  |Rx + Tx - 2C| that I first compared with only holds inside the retro-reflected beam, about
  the aperture wide; Rx_1 and Rx_2 at 70 and 140 mm lie outside it, where the far-field
  (physical-optics) path r_t + r_r applies.)
* T5: Rx_0 on top of the Tx site solves cleanly and gives the T1 level. Pass.
* First look at cross-pol (4 rays per wavelength): |S_HV|, |S_VH| relative to co-pol -60.7 dB
  (Rx_0, monostatic), -38.5 dB (Rx_1, 1 deg), -33.1 dB (Rx_2, 2 deg): about +6 dB per
  doubling of the bistatic angle, HV = VH. Pending T4 (ray density) and T4b (symmetry) before
  reading it as physics: at targetRoll 0 the trihedral's mirror plane is vertical (x-z) while
  the receivers are displaced horizontally, so cross-pol is allowed by symmetry.
* T4 ray density (4 / 8 / 16 rays per wavelength; 25 / 47 / 257 s): the monostatic cross-pol
  falls -60.7 / -75.5 / -92.8 dB, a numerical floor that keeps dropping with ray density;
  the bistatic cross-pol stays -38.2 ... -38.5 dB (1 deg) and -32.8 ... -33.1 dB (2 deg), and
  the co-pol level within 0.03 dB: converged. A repeated 4-ray solve is bit-identical. Pass;
  4 rays per wavelength suffice for anything above about -55 dB re co-pol.
* T2 polarization labels. The PEC strip (40 x 0.5 mm along H) gave |S_HH| = |S_VV| =
  -114.9 dB: physical optics ignores the polarization of a flat plate at normal incidence,
  however thin, so this test cannot tell H from V. Replaced it by a dihedral (seam along local
  x = H at targetRoll 0, `TrustDihedral`) against the trihedral, both referenced at TargetCS:
  a PEC reflection maps E -> -R(n) E (R the mirror reflection), so three bounces give
  (-1)^3 (-I) E = +E and two give a 180 deg rotation about the seam, which keeps the
  component along the seam and flips the one across it. Measured ratios dihedral / trihedral
  at Rx_0: S_HH -0.1 deg (+-6 deg over the band), S_VV -176 deg; magnitude +6.9 dB (theory
  8 pi a^4 / lambda^2 against 4 pi a^4 / (3 lambda^2): +7.8 dB). "Horizontal" is along y and
  "Vertical" along z, as assumed. Pass.
* T4b symmetry (`t4b_trihedral_roll90_1deg`): at targetRoll 90 deg the bistatic plane (x-y)
  is the trihedral's mirror plane, and the bistatic cross-pol drops from -38 / -33 dB to
  -60 / -59 dB, the 4-ray numerical floor; the co-pol level is unchanged. Pass. So the
  cross-pol at targetRoll 0 is genuine bistatic depolarization of the trihedral (no mirror
  symmetry about the bistatic plane), growing about 6 dB per doubling of the bistatic angle
  (amplitude proportional to beta): the first finding of the action point, to be swept in M4.
* Trust checks done (all pass); stopped for validation before the production sweep. The
  sandbox designs: `MultiRxTrihedral` (reset to targetRoll 0), `TrustStrip`, `TrustDihedral`;
  exports `t1_*`, `t2_*`, `t4_*`, `t4b_*` in the sandbox `exports/`.
* Validation (user): production approved with a dihedral at 22.5 deg added; sphere 40 mm across
  (radius 20 mm) like the other targets; per-receiver analysis mode approved; horizontal
  displacement only, vertical noted as an extension action point.
* Mini sweep (`MultiRxSphere`, bistaticAngle 1 and 2 deg x rxAim 0 and 1; 4 solves in 1.1 min)
  checked the sweep -> export -> pooling chain: rxAim changes |S| by at most 0.05 dB at 4 deg
  (the 60 deg beams' roll-off), and Rx_1 at 2 deg equals Rx_2 at 1 deg to 7e-9 (same position).
  It also exposed the sphere's faceting (T6 above), so the sphere runs at mesh slider 9.
  (Slider 10 is rejected by `InitialMeshSettings`.)
* Production launched (`run_production.sh`).
* Sphere done (64 solves, 23.8 min; 49 distinct bistatic angles 0 ... 16 deg after pooling;
  equal-angle points agree to 3e-8). Co-pol invariant: |S_VV / S_HH| within -0.13 ... +0.11 dB
  at 79 GHz (band -0.61 ... +0.66 dB), co-pol phase difference -4.1 ... +3.2 deg up to 16 deg;
  H <= 0.006, mean alpha <= 0.5 deg. Cross-pol at the numerical floor (co/cross median 98 dB),
  so every cross-pol quantity (co/cross ratio and phase, HV/VH, asymmetry, k3, k4) is noise and
  its phase scatters at random: an erratic cross-pol phase means nothing while the cross-pol
  sits at the floor. (Refined 2026-10-24: at a fixed angle the floor's phase is smooth in
  frequency, not random; it is meaningless because it changes with ray density.) The span falls 0.4 dB by 16 deg with aimed receivers (the longer receive
  path, R / cos 16 deg: 0.34 dB) and 1.2 dB on boresight (plus the 60 deg beams' roll-off);
  rxAim changes no polarimetric quantity, so in this principal plane the parametric beams'
  off-axis polarization definition adds no squint.
* Trihedral done (128 solves, 62 min). Departure angles from `summary_departures_aim1.csv`
  (receivers aimed; first bistatic angle, at 79 GHz / anywhere in the band):

  | Descriptor | Level | Roll 0 (mirror plane vertical) | Roll 90 (mirror plane = bistatic plane) |
  | --- | --- | --- | --- |
  | span (retro lobe) | -3 dB / -10 dB | 3.5 / 5.75 deg | 3.5 / 6.0 deg |
  | co-pol ratio abs(VV/HH) | 0.5 dB | 7.5 (6.75) deg | 6.25 (2.25) deg |
  | co-pol phase difference | 5 deg / 10 deg | 4.5 / 5.75 deg | 5.75 / 6.5 deg |
  | co/cross per Tx pol | below 30 dB / 20 dB | 2.25 ... 2.5 / 7.0 ... 7.5 deg | 8.5 deg / never (floor) (superseded below) |
  | entropy H | 0.1 / 0.3 | 7.5 / 7.75 deg | 7.5 / 8.5 deg |
  | mean alpha | 10 deg / 20 deg | 6.5 / 7.25 deg | 6.75 / 7.5 deg |

  The first null of the bistatic pattern lies at 7.8 deg (roll 0) and 8.5 deg (roll 90; the
  triangular aperture is wider in that plane); lambda / D = 7.8 deg gives an effective
  aperture D = 28 mm. Inside the main lobe (below about 6 deg) the co-pol descriptors and
  Cloude-Pottier stay put; they all break down together at the null, where the triple bounce
  vanishes and edge and single/double-bounce contributions compete (H = 0.48, mean alpha
  43 deg at roll 0).
* Trihedral cross-pol: at roll 0 it rises out of the floor at once and is the first
  descriptor to move: co/cross falls below 30 dB by 2.25 deg, i.e. inside the -3 dB retro lobe,
  with a deterministic co/cross phase (about +135 deg for both Tx polarizations, flat over the
  band and over 0.25 ... 6 deg; corrected 2026-10-24: +135 deg is the old convention
  arg(S_HH S_VH*) at 79 GHz only, and the phase is not flat but linear over the band) and HV = VH within 1 dB up to about 4 deg (reciprocity holds; the
  asymmetry and k4 grow from there). At roll 90 the cross-pol stays at the numerical floor
  (about -145 dB absolute) at every angle; its co/cross and HV/VH "departures" near 7.5 ... 9
  deg are the co-pol null falling onto that floor, not a cross-pol: the 50 dB co/cross guard of
  the departure table is relative to co-pol and cannot tell the two apart near a null.
* So for a trihedral the order of breakdown is: cross-pol (if the bistatic plane is not a mirror
  plane) from about 1 deg, co-pol phase difference at about 4.5 deg, then co-pol ratio,
  entropy and alpha together at the edge of the main lobe (6.5 ... 7.5 deg).
* Dihedral done (320 solves, 151 min; 10 exports). The departure table is now relative to each
  target's own monostatic value. Its earlier absolute levels (co/cross below 30 dB, mean alpha
  above 10 deg) assumed a co-pol-pure, odd-bounce target and fired at the first step for every
  dihedral (mean alpha 90 deg; co/cross about 0 dB at roll 22.5 and 67.5 deg, -38.6 dB at 45 deg).
  Kinds now: 'change' from the monostatic value; 'floor change' for co/cross (clipped to
  +-50 dB, so a cross-pol rising out of the floor counts from there; 20 dB = the former
  "below 30 dB"); 'absolute' / 'above' only for HV/VH and the asymmetry, which reciprocity
  fixes monostatically. A descriptor whose monostatic value is a floor or a reciprocity zero is
  'undefined' (and is left out of `summary_change_<family>_aim<a>.png` per design).
* The floor masks are now absolute: each cross-pol (co-pol) channel must lie within 50 dB of
  the target's monostatic span. The floor is fixed relative to that span, not to the local
  co-pol: about -60 dB for the trihedral and -85 ... -100 dB for the sphere and the dihedral,
  flat over 0 ... 16 deg (trihedral roll 90: HV -63.7 ... -59.4 dB, VH -63.7 ... -61.6 dB).
  This removes the trihedral roll 90 caveat above: its late co/cross and HV/VH "departures" are
  gone (undefined throughout), as they were a co-pol null falling onto a fixed floor.
* Dihedral checks. Monostatic span -75.3 dB at every roll (spread 0.02 dB: roll-invariant),
  35.0 dB above the sphere (theory 8 pi a^2 b^2 / lambda^2 against pi r^2: 35.5 dB) and 8.0 dB
  above the trihedral (7.8 dB). Pauli fractions follow the roll: k2 = 0 dB at 0 and 90 deg;
  k2 = k3 = -3.0 dB at 22.5 and 67.5 deg (co/cross +-0.14 dB, the expected mixed response);
  k3 = 0 dB at 45 deg. k1 = -38.6 dB at every roll is the dihedral's own HH/VV imbalance
  (|VV/HH| -0.19 dB, phase 179.5 deg, from edge diffraction); the trace is roll-invariant, so at
  45 deg it is the residual co-pol (co/cross -38.6 dB, HH = VV in phase). 67.5 deg mirrors
  22.5 deg monostatically (co/cross phase H -179.7 / V -0.3 deg against -0.3 / -179.7 deg).
  Mean alpha 89.6 deg at every roll.
* Dihedral lobes. First null 5.75 deg at roll 0 (seam along H, in the bistatic plane: the
  dihedral is a mirror along its seam, lambda / L = 5.44 deg for L = 40 mm) and 4.0 deg at roll
  90 (seam across the bistatic plane: retro-reflecting, lambda / (2 x 40 mm cos 45 deg) = 3.85
  deg); 6.0, 5.5 and 4.25 deg at 22.5, 45 and 67.5 deg. So bistatically 67.5 deg does not mirror
  22.5 deg: the lobe width follows the seam's orientation to the bistatic plane.
* rxAim (all dihedrals, as for the sphere and the trihedral) changes the span only (up to
  +0.85 dB at 16 deg: the 60 deg beams' gain) and no polarimetric quantity outside nulls.
* Departures, all targets (receivers aimed; first bistatic angle at 79 GHz, in units of the
  first null beta_0 in brackets; '-' never up to 16 deg, 'undef.' monostatic value a floor):

  | Descriptor (level) | D 0 | D 22.5 | D 45 | D 67.5 | D 90 | T 0 | T 90 |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | first null beta_0 | 5.75 | 6.0 | 5.5 | 4.25 | 4.0 | 7.8 | 8.5 |
  | span (-3 dB) | 2.75 (0.48) | 2.5 (0.42) | 2.25 (0.41) | 2.0 (0.47) | 2.0 (0.50) | 3.5 (0.45) | 3.5 (0.41) |
  | co-pol ratio (0.5 dB) | 5.25 (0.91) | 5.5 (0.92) | 1.0 (0.18) | 2.75 (0.65) | 2.75 (0.69) | 7.5 (0.96) | 6.25 (0.74) |
  | co-pol phase (5 deg) | 4.75 (0.83) | 4.75 (0.79) | 1.25 (0.23) | 4.0 (0.94) | 4.0 (1.0) | 4.5 (0.58) | 5.75 (0.68) |
  | co/cross, Tx H / V (3 dB) | - | 6.25 / 6.0 (1.0) | 2.25 / 1.5 (0.27 ... 0.41) | 4.25 (1.0) | - | 0.5 (0.06) | - |
  | co/cross, Tx H / V (20 dB) | - | - | 5.25 (0.95) | - | - | 2.5 / 2.25 (0.29 ... 0.32) | - |
  | co/cross phase, Tx H / V (5 deg) | undef. | 4.75 / 5.25 (0.79 ... 0.88) | 1.25 / 2.0 (0.23 ... 0.36) | 3.5 / 4.0 (0.82 ... 0.94) | undef. | undef. | undef. |
  | HV/VH (1 dB) | - | 5.75 (0.96) | 5.25 (0.95) | 4.25 (1.0) | - | 0.5 (0.06) | - |
  | asymmetry (-30 dB) | - | 4.5 (0.75) | 3.5 (0.64) | 3.5 (0.82) | - | 7.5 (0.96) | - |
  | entropy H (0.1) | - | 6.0 (1.0) | 5.25 (0.95) | 4.25 (1.0) | 3.75 (0.94) | 7.5 (0.96) | 7.5 (0.88) |
  | mean alpha (10 deg) | 5.75 (1.0) | 13.0 (2.2) | 5.25 (0.95) | 4.25 (1.0) | 3.75 (0.94) | 6.5 (0.83) | 6.75 (0.79) |

  The sphere departs nowhere (co-pol ratio 0.5 dB only somewhere in the band at 14 deg).
  The trihedral roll 0 HV/VH crosses 1 dB at 0.5 deg, where its cross-pol is only about 16 dB
  above the trihedral's floor; the 3 dB level is crossed at 5.25 deg.
* Two regimes. (1) Descriptors built from channels that are strong monostatically hold to
  0.6 ... 1.0 beta_0 (0.5 dB, 5 deg, 3 dB co/cross with mixed co and cross at 22.5 / 67.5 deg,
  entropy 0.1, alpha 10 deg), well beyond the span's -3 dB point at 0.41 ... 0.50 beta_0; they
  go together at the first null. (2) Channels that are zero or weak monostatically move first:
  the trihedral's cross-pol at roll 0 (co/cross 3 dB down at 0.06 beta_0, 20 dB down at 0.3
  beta_0) and the 45 deg dihedral's residual co-pol (0.5 dB / 5 deg at 0.18 ... 0.23 beta_0).
  Where the bistatic plane is a mirror plane of the target (trihedral 90, dihedral 0 and 90,
  sphere), cross-pol never appears.
* Order of the emerging channels in beta (fit over 0.25 ... 1.5 deg, 79 GHz): the trihedral's
  cross-pol at roll 0 grows 5.1 dB per doubling of beta (-49 / -44 / -39 / -34 dB re the monostatic span
  at 0.25 / 0.5 / 1 / 2 deg): field proportional to beta. The rolled dihedrals' emerging
  channels grow 12 ... 14 dB per doubling: k2 at 45 deg (-82 / -68 / -55 dB at 0.5 / 1 / 2 deg)
  and k4 at 22.5 and 45 deg: field proportional to beta^2. The trihedral's k4 stays at its floor
  (-60 dB) below 2 deg: its first-order cross-pol is symmetric (HV = VH).
* T7 (why the order differs): turning the scene by 180 deg about the Tx-target line maps the
  receiver at +beta to -beta and (H, V) to (-H, -V), which leaves every S_pq unchanged. A target
  with that C2 symmetry (any dihedral, plate or sphere at boresight) therefore has
  S(-beta) = S(beta): even in beta, with no linear term. A mirror plane x-z instead gives
  S(-beta) = M S(beta) M, M = diag(-1, 1): co-pol even, cross-pol odd, so a linear cross-pol
  is allowed. The C3 trihedral at roll 0 has that mirror plane but no C2. Solved with the
  receivers mirrored to -y (bistaticAngle -2 deg, rxAim 1; exports `mirror_*`):
  dihedral 22.5 deg: S(-beta) / S(beta) = 0.000 dB, 0.00 deg on every channel (bit-identical
  ray tracing); dihedral 45 deg: within 0.03 dB / 0.2 deg (the -45 dB weak co-pol channels);
  trihedral 0 deg: co-pol within 0.03 dB / 0.15 deg, cross-pol sign-flipped (177.7 ... 179.9
  deg), residuals S(-beta) -+ S(beta) at -50 ... -68 dB re peak, the trihedral's floor. Pass.
  (These solves left the designs' nominal values at bistaticAngle -2 deg, rxAim 1.)
* Earlier note findings, explained: the erratic cross-pol phase (progress meeting) belongs to
  cross-pol at the numerical floor (sphere; any target whose mirror plane is the bistatic
  plane); a genuine bistatic cross-pol has a deterministic phase (trihedral roll 0: about +135
  deg, flat over the band; rolled dihedrals: 0 / 180 deg). (Corrected 2026-10-24: the
  trihedral's phase is linear over the band, and a floor's phase is smooth too; what tells them
  apart is convergence with ray density.) Cloude-Pottier's stability holds
  (H and alpha move only at 0.8 ... 1.0 beta_0) but is partly blindness: at roll 0 the
  trihedral's co/cross had fallen 20 dB by 2.5 deg while H and alpha stayed put to 6.5 deg.
* Scaling check (user: "Run it"; `run_scaling.sh`, 60 solves, 92 min, log `scaling.log`;
  analysis `scaling.py` -> `scaling/scaling.png`, `scaling.csv`). Trihedral at roll 0:
  40 mm at 1 m (the 4 m run exists), 80 mm at 4 m and 16 m, 160 mm at 4 m; dihedral at roll 90
  (retro plane): 160 mm at 4 m. bistaticAngle grids scaled with the expected first null; rxAim 1.
  Far-field distance 2 (sqrt 2 a)^2 / lambda: 1.7 m (40 mm), 6.7 m (80 mm), 27 m (160 mm).
* First null: beta_0 = 8.03 / 7.79 deg (40 mm at 1 / 4 m), 3.95 / 3.97 deg (80 mm at 4 / 16 m),
  1.99 deg (160 mm): proportional to 1 / a within 2 %, independent of range, also inside the
  target's far-field distance. The span collapses onto one curve against beta / beta_0 out to
  the first null (beyond it the 160 mm trihedral departs).
* Strong-channel departures in beta_0 units hold across sizes and ranges: span -3 dB at
  0.44 ... 0.50, co-pol ratio 0.5 dB at 0.62 ... 1.0, co-pol phase 5 deg at 0.56 ... 0.70 (160 mm:
  3.0), entropy 0.1 at 0.96 ... 1.0, mean alpha 10 deg at 0.83 ... 1.0. So the lobe rule stands:
  the co-pol signature is preserved while beta < about 0.5 beta_0, i.e. d < 0.5 R lambda / D.
* Trihedral cross-pol (re monostatic span, 79 GHz), the hypothesis "a function of beta / beta_0"
  is falsified. It is a function of beta alone at a given size: 40 mm at 1 m against 4 m
  -44.3 / -44.0 dB at 0.5 deg, -38.6 / -38.8 dB at 1 deg, -34.0 / -34.3 dB at 2 deg; 80 mm at 4 m
  against 16 m within 0.5 dB. It falls with size: at 1 deg -38.8 / -42.7 / -57.2 dB (HV) for
  40 / 80 / 160 mm, peak inside the lobe -33.2 / -41.3 / about -52 dB, about -9 dB per doubling.
  Its shape against beta / beta_0 is common (rising about 4 ... 5 dB per doubling, peak near
  0.4 beta_0); only the level moves.
* T8 edge diffraction (`TrustTrihedralEdges`, 40 mm at roll 0, bistaticAngle 1 deg; exports
  `t8_*`): HV re span at 1 / 2 deg -49.4 / -44.8 dB with physical optics only, -41.7 / -37.2 dB
  with PTD, -38.8 / -34.3 dB with PTD + UTD rays (production). The bistatic cross-pol exists in
  physical optics, but edge diffraction adds about 10 dB, so its level is model-dependent
  (consistent with its fall with size, edges against aperture). The monostatic cross-pol is
  -142 dB without UTD rays and -63.7 dB with them: the trihedral's "floor" is UTD-ray noise
  (which is why it fell with ray density in T4). Co-pol within 0.3 dB in all three models.
* Dihedral 160 mm at 4 m (27 m far-field distance): no lobe and no null in 0 ... 6 deg. The span
  rises +3.1 / +5.4 dB at 0.5 / 1 deg and stays within +-3 dB of monostatic up to 2.5 deg (a
  Fresnel ripple over the geometric footprint of the retro-reflected beam, half-angle about
  atan(226 mm / 4 m) = 3.2 deg), then falls (-9 dB at 4 deg). Its polarimetric ratios stay put:
  |VV/HH| within -0.23 ... +0.43 dB and the co-pol phase within 180 +- 4.4 deg up to 4 deg. Its
  cross-pol stays at the floor (symmetric orientation). Near-field: the power varies across an
  array, the polarization ratios do not.
* Correction of my earlier extrapolation to the user: "-28 dB for a 15 cm calibration trihedral
  at 3 m" assumed the beta / beta_0 scaling and is wrong. Interpolated between 80 and 160 mm, a
  15 cm trihedral's bistatic cross-pol at 1 deg is about -55 dB re span (PTD + UTD; about
  10 dB lower in physical optics), well below a real antenna's cross-pol.

### 2026-10-21

* Merged: `main`'s uncommitted simulator changes became six topical commits (linked discrete
  sweeps, mapped variables, builder presets, robust exports, x boresight, co/cross signatures),
  the branch was rebased onto them (same final tree) and merged; the worktree was removed.
  `SANDBOX_DIR` now follows `pomaar.config.INTERIM_DATA_DIR` instead of a home-relative path.
* Receiver aiming closed as a reference (`investigation/aim_check.py`, outputs and argument in
  `rx_aim/`). Per receiver, angle and frequency, S_aim = g S_boresight:
  * The residual of that fit is at most -48 dB inside the lobe (span within 20 dB of
    monostatic). Its median there is -79 ... -103 dB per design, and it stays at most -32 dB
    even in nulls.
  * |g|^2 = 3 (beta / 30 deg)^2 dB (the 60 deg beam) to within 0.1 dB, and |arg g| < 0.8 deg.
  * Since every descriptor but the span is invariant under S -> g S, aiming changes the span
    only.
  * Scope: identical cross-pol-free H/V beams in a principal plane. With real patterns, aiming
    changes the one-way receive polarization matrix P_Rx(beta): the E/H-plane gain ratio, the
    pattern cross-pol and, off the principal planes, the basis rotation. That comes from the
    far-field patterns, not from SBR+.
* Signature conventions switched to the radar ones, under new keys so that no old name changes
  meaning. Each new quantity is the negative of the old one (in dB or deg), and the departure
  angles are unchanged (all 192 rows of `summary_departures_aim1.csv` checked):

  | Old key (CSVs before 2026-10-21) | Old quantity | New key | New quantity |
  | --- | --- | --- | --- |
  | `copol_ratio_db` | \|S_VV/S_HH\| | `zdr_db` | Z_DR = \|S_HH\|²/\|S_VV\|² |
  | `copol_phase_difference_deg` | ∠(S_VV S_HH*) | `hh_vv_phase_deg` | ∠(S_HH S_VV*) |
  | `co_minus_cross_txh_db` | \|S_HH/S_VH\| | `ldr_txh_db` | LDR_H = \|S_VH\|²/\|S_HH\|² |
  | `co_minus_cross_txv_db` | \|S_VV/S_HV\| | `ldr_txv_db` | LDR_V = \|S_HV\|²/\|S_VV\|² |
  | `co_cross_phase_txh_deg` | ∠(S_HH S_VH*) | `cross_co_phase_txh_deg` | ∠(S_VH S_HH*) |
  | `co_cross_phase_txv_deg` | ∠(S_VV S_HV*) | `cross_co_phase_txv_deg` | ∠(S_HV S_VV*) |

  The "co/cross" levels quoted in this log (e.g. "co/cross falls below 30 dB") read as
  LDR above -30 dB in the new convention. The analysis figures are now split by group (powers,
  copolar, crosspolar, pauli, cloude_pottier) with mathtext labels.
* Labels (user: "'Tx H: |S_HH| - |S_VV|' is mathematically unclear"): power ratios are written
  as ratios on a dB axis (`LDR_H = |S_VH|^2/|S_HH|^2`, axis in dB), never as differences of
  magnitudes; phases as the argument of a product (`arg(S_VH S_HH*)`). The shared labels are
  `bistaticity_analysis.LABELS`, used by `summarize.py` and `scaling.py`.
* M7 validation (`investigation/validation/README.md`, all on a copy of the project). Bugs fixed:
  PyAEDT creating an empty design for a missing name (builder, `sbr_simulator`, `sbr_results`);
  `--target none` unsolvable (SBR+ needs a model object; now a vacuum placeholder); failed solves
  exiting with 0 (now exit code 3); `--overwrite` hanging on the 'asol_priv in use' modal (now
  save, close and reopen); invalid `--map-variables` caught only after building; `--scene-only`
  building an x-style TargetCS in a z-boresight design.
* The z-boresight quirk, measured (S3/S4: dihedral at roll 45 deg, 10 deg in azimuth and 5 deg
  in elevation, x and z scenes). In a z scene `targetAzimuth`/`targetElevation` are spherical
  angles (about and off boresight; here 26.74 / 11.17 deg), while the processor's angles run
  towards H and towards V (10 / 5 deg). Steered at the right direction, both scenes give the same
  cross-pol phase arg(S_HV S_VH*): +3.7 / +3.6 deg (from -17.9 deg at boresight; the 21.7 deg
  removed equals the prediction from the virtual positions). Steered with the scene's spherical
  angles: +44.4 deg. It only mattered off boresight, in z scenes, when the scene's angles were
  passed to the processor. Exports now carry `target_steering_deg` (processor angles) and
  `boresight`, and `bistaticity_analysis` steers towards it by default. All bistaticity scenes
  here are x-boresight with the target on boresight, so no earlier result is affected.
* Summary bug: `summarize.figure` paired series with colours by `zip`, so it silently dropped
  every series after the eighth. The earlier `summary_signatures_aim1.png` and
  `summary_change_aim1.png` (13 designs) lacked the trihedral rolls 0 / 90 and TrihedralS160;
  the CSVs and the departure tables (and so every number in this log) were complete. The
  figures are now one per target family (`summary_*_<family>_aim<a>.png`, at most eight series,
  more raise an error).
* M8 extension run (`run_extension.sh`, log `extension.log`, 74 min; 20 exports):
  `MultiRxTrihedral` at targetRoll 45 (slant), `MultiRxPlate` (40 mm square PEC plate, normal
  to boresight), and the split-pair scenes `PairsSphere`, `PairsTrihedral` (roll 0, 90),
  `PairsDihedral` (roll 0, 22.5), pairAngle 0.5 ... 6 deg, gapRatio 0 and 2, all at 4 m.
* Slant trihedral (roll 45; departures from `summary_departures_aim1.csv`, 79 GHz / band):
  * span: as roll 0 (-3 dB at 3.5 deg, -10 dB at 5.75 deg): the lobe follows the aperture,
    whose width hardly changes with roll;
  * cross-pol: rises at once like roll 0, first order (LDR_H -51 / -45 / -39 / -33 dB at
    0.25 / 0.5 / 1 / 2 deg; roll 0: -50 / -43 / -36.5 / -31 dB), 1 ... 2.5 dB lower;
  * co-pol: unlike roll 0, it changes to first order too (arg(S_HH S_VV*) -0.5 / -1.2 / -2.6 deg
    at 0.5 / 1 / 2 deg, Z_DR -0.1 / -0.2 / -0.3 dB), so Z_DR departs 0.5 dB at 3.5 deg against
    7.5 deg at roll 0, and the co-pol phase 5 deg at 3.5 against 4.5 deg;
  * entropy and mean alpha: as roll 0 (0.1 at 7.0 deg; 10 deg at 6.25 deg).
* That completes the symmetry picture (T7): with the displacement along H, a mirror plane x-z
  (perpendicular to the displacement) makes co-pol even in beta (second order) and cross-pol
  odd; a mirror plane x-y (the bistatic plane) forbids cross-pol but leaves the co-pol free to
  change linearly; neither (slant trihedral) leaves both linear. Trihedral roll 0 (x-z):
  co-pol phase -0.2 / -0.8 / -3.9 deg at 1 / 2 / 4 deg (quadratic); roll 90 (x-y): +0.3 / +0.6 /
  +1.1 / +1.6 deg at 0.25 / 0.5 / 1 / 2 deg and Z_DR +0.1 / +0.2 / +0.3 dB at 0.5 / 1 / 2 deg
  (linear), no cross-pol; roll 45: both linear.
* Flat plate (40 mm, normal to boresight): monostatic span -77.7 dB, 2.4 dB below the 40 mm
  dihedral (PO: 4 pi A^2 / lambda^2 against 8 pi a^2 b^2 / lambda^2, 3 dB). Its bistatic
  pattern is the plate's sinc^2 with the first null near lambda / L = 5.4 deg (span -10 dB at
  4.25 deg, -20 dB by 6 deg), like the dihedral at roll 0, whose aperture along H is the same
  40 mm. Departures likewise close to that dihedral: Z_DR 0.5 dB at 4.75 deg (dihedral 5.25),
  co-pol phase 5 deg at 4.0 deg (4.75), mean alpha 10 deg at 5.25 deg (5.75). Its cross-pol is
  exactly zero in physical optics (LDR -150 ... -360 dB, never defined): a flat PEC plate in a
  principal plane does not depolarize.
* Split-pair Z_DR (`pairs.py` -> `pairs/pairs_zdr.png`, `pairs_zdr.csv`). The
  supervisor's layout TxV -[b]- RxV -[g]- TxH -[b]- RxH: VV and HH each from a pair of base b
  (pair bistatic angle beta), the pairs g apart, so their bisectors lie (b + g) / R apart
  (beta and 3 beta for g = 0 and 2 b). Reference: a co-located dual-polarized pair of the same
  base (bisector on boresight). Split minus reference, at 79 GHz (band max |.| in brackets):

  | Target | g | 1 deg | 2 deg | 3 deg | 6 deg |
  | --- | --- | --- | --- | --- | --- |
  | sphere | 0 / 2b | 0.01 / 0.00 dB | 0.12 / 0.10 dB | 0.13 / -0.02 dB | -0.04 / -0.09 dB (0.12) |
  | trihedral roll 0 | 0 / 2b | -0.03 / 0.00 dB | -0.05 / -0.04 dB | -0.02 / -0.03 dB | 0.72 (0.94) / 0.18 (1.46) dB |
  | trihedral roll 90 | 0 / 2b | 0.03 / 0.11 dB | -0.05 / -0.18 dB | -0.10 (0.54) / -1.26 (1.40) dB | -2.54 (6.25) / -7.55 (12.6) dB |
  | dihedral roll 0 | 0 / 2b | 0.05 / -0.01 dB | 0.01 / 2.0 (8.6) dB | -0.02 / 0.22 dB | 3.0 (12.2) / 3.7 (4.0) dB |
  | dihedral roll 22.5 | 0 / 2b | 0.05 / -0.14 dB | -0.05 / 4.5 (15.5) dB | -0.17 / -0.34 (1.39) dB | 0.76 (47.9) / 3.2 (6.3) dB |

  The HH and VV pair angles match exactly (mirror images). The co-pol phase difference
  arg(S_HH S_VV*) changes by at most 0.7 deg up to 2 deg with g = 0 and 2.4 deg at 1 deg with
  g = 2b; the dihedral at 22.5 deg reaches -47 deg at 2 deg with g = 2b and -61 deg at 6 deg
  with g = 0.
* Reading of the split pairs. The VV pair is the HH pair mirrored about x-z, so the split Z_DR
  compares the target seen from two mirrored aspects, (b + g) / R apart. For a target with a
  mirror plane x-z or C2 symmetry about the line of sight (sphere, trihedral roll 0, dihedrals)
  the two aspects are equivalent and the split error reduces to the change of Z_DR with the
  bisector's offset, which is even (second order). The trihedral at roll 90 has neither and gives
  the largest, steadily growing error: -1.3 dB at 3 deg and -7.6 dB at 6 deg with g = 2b,
  against a reference pair that changes by only 0.3 / 1.4 dB. The dihedrals' large values with
  g = 2b at 2 and 4 deg are isolated spikes: there the pairs' bisectors sit +-3 / +-6 deg off
  boresight, near the 40 mm aperture's monostatic aspect nulls (lambda / (2 L) = 2.7 deg and its
  multiples), where both co-pol channels are small and their ratio is ill-conditioned. Below
  2 deg with g = 0 (bisectors under 2 deg apart) every target stays within 0.2 dB (band) and
  0.7 deg of the co-located pair: at automotive geometry (b of a few cm at a few metres,
  beta below 1 deg) the split costs nothing with these ideal beams. The order of growth below
  2 deg is not resolvable: the differences there sit at the ~0.05 dB level of the band spread.
* Results integrated into the main tree (user: "the results are solid"; then "only an actual
  report ... will sit [in reports]" and "scrap the whole sandbox folder"). The sandbox folder is
  gone:
  * the AEDT project moved to
    `~/Projects/AEDT/bistaticity-parametric-beams/bistaticity_parametric_beams.aedt` (renamed
    with its results folder; the scripts' default, `investigation/project.py`, and `--project`
    for another);
  * the exports and run logs moved to `data/interim/bistaticity/parametric_beams/`, and the
    results to `data/processed/bistaticity/parametric_beams/` (untracked, regenerable);
  * this log, the README and the receiver-aiming argument moved next to the scripts in
    `src/pomaar/simulator/investigation/`;
  * the validation copy moved to
    `~/Projects/AEDT/mimo-polarimetry-validation/mimo_polarimetry_validation.aedt`, with its
    port tables, and its logs and exports to `data/interim/validation/`.

  The `parametric_beams` layer marks the idealized antennas, next to the full-wave cluster
  results in `data/*/bistaticity/`. Every regenerated value matched the log (e.g. T1 +0.25 dB,
  the scaling's beta_0 and HV levels).
  A re-export from the moved project (`MultiRxPlate`) is bit-identical to the original.
  The pooling's consistency check reports 32 deg between equal-angle points of the plate, but
  only in its cross-pol phases: those channels are exactly zero in physical optics (−256 …
  −443 dB, double-precision rounding), so their phases are arbitrary. The plate's co-pol, like
  every other export's, agrees to 1e-7.

### 2026-10-24

* Presentation figures (`investigation/slides.py`, one message per figure) written to
  `reports/bistaticity-analysis/parametric_beams/`; the horn-cluster exports re-analysed with
  the current `bistaticity_analysis` (radar conventions, figures by group) into
  `reports/bistaticity-analysis/horn_cluster_*`.
* Correction, trihedral cross–co phase. Drawn against frequency, arg(S_VH S_HH*) of the
  trihedral at roll 0 is not flat but linear over the band. A linear fit over 74–84 GHz gives a
  two-way path difference of 46.8 … 47.7 mm (40 mm edges, 0.5 … 4 deg), 46.0 … 46.2 mm (roll
  45), 91.9 … 92.3 mm (80 mm) and 183.5 mm (160 mm, 1 deg), with 6 … 12 deg rms about the line.
  That is 2a/√3 (46.2 / 92.4 / 184.8 mm), twice the depth of the corner below the aperture
  plane: the bistatic cross-pol comes from the aperture rim, the co-pol from the corner. The
  positive slope (shorter path in AEDT's e^{+jωt}) puts the rim in front of the corner, and the
  edge origin fits T8 (edge diffraction adds 10 dB). At 79 GHz the phase is −124 … −134 deg;
  the "+135 deg" of the earlier entries was the old convention arg(S_HH S_VH*) at that frequency.
  The rolled dihedral (22.5 deg) has no such delay (0.1 mm, 0.8 … 0.9 deg rms, −179.5 deg): its
  cross-pol is the same double bounce as its co-pol.
* Refinement, floor phase. At a fixed bistatic angle the floor's phase is a smooth function of
  frequency (deterministic ray tracing), not random; the "erratic" look came from sweeps and
  band spreads. What separates a floor from a genuine cross-pol is convergence with ray density
  (T4 exports, trihedral roll 0, 79 GHz): monostatic LDR_H −60.7 / −74.2 / −103.3 dB at 4 / 8 / 16
  rays/λ with unrelated phases; 1 deg −36.5 / −35.8 / −36.0 dB and 2 deg −31.0 / −30.6 / −30.7 dB,
  phases within about 8 deg over the band (`07_cross_pol_genuine_vs_floor.png`). The horn
  cluster's erratic cross-pol is a third case: genuine antenna cross-pol of different element
  pairs, −33 … −73 dB re co-pol, changing with spacing.

### 2026-10-28

* `single_pol_power_only.py` (architecture A, single-pol power-only, against B, dual-pol
  coherent): Pauli fractions now of the span (f = |x|^2 / span; f^ over the span from the four
  H/V powers, eight measured powers), both scenes on one target list (`TARGETS`) with
  per-design figures. The three-receiver worst case exceeds 0 dB (+28 … +48 dB) because the
  pairs at 0, beta, 2 beta see target powers up to 30 … 39 dB apart (aperture-lobe nulls), not
  because the estimated total cancels (earlier explanation, corrected).
* Split-pair runs for the targets without a split-pair scene (`run_split_pairs.sh`, log
  `logs/split_pairs.log`). Pair B of the trihedral (no symmetry mapping +y onto -y at roll 45
  or 90 deg) comes from the receiver ladder at roll -psi (270, 315 deg) mirrored in x-z, the
  trihedral's mirror plane at roll 0. The project would not open at first: a stale
  `.aedt.lock` of a dead AEDT process (21 Oct), moved aside. The run was aborted by the user
  during the trihedral ladders. Done and exported: `PairsDihedral_gapRatio0_targetRoll45`,
  `PairsTrihedral_gapRatio0_targetRoll45`, `PairsPlate_gapRatio0` (new design `PairsPlate`),
  `MultiRxPlate_rxAim0`. Solved, not exported: `MultiRxTrihedral` setup `Rx0Roll45Sweep`
  (rxAim 0, roll 45 deg). Failed or aborted (AEDT reported an error after 0.2 … 0.4 min):
  `Rx0Roll270Sweep`, `Rx0Roll315Sweep`. `sweep.py` printed "Solved" anyway; it now raises
  `SolveFailedError` when `analyze_setup` fails.

### 2026-10-29

* Finished the split-pair run (the rest of `run_split_pairs.sh`): `Rx0Roll270Sweep` and
  `Rx0Roll315Sweep` solved (12.3, 11.6 min); exported `MultiRxTrihedral_rxAim0_targetRoll{45,270,315}`.
* Scope correction (Martin): the split pairs answer only the supervisor's question, Z_DR from
  TxV–RxV and TxH–RxH pairs of the same base with shifted bisectors (`pairs.py`); the
  architecture comparison (A power-only against B coherent) belongs to the three-receiver
  scenes, which measure the bistatic angle without assuming a layout. Removed the split-pair
  part of `single_pol_power_only.py` and its results. It had borrowed the full S of pairs A
  and B from the receiver ladder (`MultiRx*_rxAim0`, by reciprocity and the target's symmetry
  or mirror image), i.e. dual-pol data from the other scene; its numbers of 28/29 Oct
  (split-pair "worst assignment" errors, the mirrored-ladder check) are withdrawn. The ladders
  at roll 270 and 315 deg served only that reconstruction.
* `pairs.py` rerun with the new split-pair scenes (dihedral and trihedral at roll 45, plate;
  gap 0 only). Split minus co-located Z_DR for beta <= 2 deg at g = 0 (79 GHz, band in
  brackets): plate 0.01 dB / 0.84 deg, trihedral roll 45 0.08 (0.18) dB / 0.51 deg, the earlier
  targets unchanged (<= 0.12 dB, <= 0.68 deg). Dihedral roll 45: 0.94 (2.13) dB / 8.8 deg, but
  its co-pol (16 … 44 dB below its cross-pol over 0.5 … 6 deg) is not reciprocal in SBR+ (VV
  of pair A against the ladder's VV, transposed: up to 13 dB apart), so its Z_DR is a ratio of
  two channels the solver does not resolve.
* Also dropped from `single_pol_power_only.py`: the three-receiver "worst assignment" column
  (A's powers spread over the receivers at 0, beta, 2 beta). It modelled a non-co-located layout
  inside the bistatic-angle study; A's powers now always come from the same Tx-Rx pair as B's.
  Deleted the exports `MultiRxTrihedral_rxAim0_targetRoll{270,315}` (only for the dropped
  reconstruction; the solutions stay in the AEDT project as `Rx0Roll270Sweep`/`Rx0Roll315Sweep`).

## Distillation (for the task note)

### Findings: multiple receivers, boresight canonical targets at 4 m (2026-10-20)

Scene: ideal parametric beams (60 deg, dual-pol: H and V share each site), one Tx and three Rx
sites at bistatic angles 0, beta, 2 beta (beta = 0.25 ... 8 deg, so 0 ... 16 deg in all),
horizontal displacement only; sphere, trihedral (roll 0, 90 deg) and dihedral (roll 0, 22.5,
45, 67.5, 90 deg), all 40 mm across. 512 SBR+ solves, sandbox branch `sandbox/bistaticity`.

* Strong-channel descriptors are lobe-limited, not displacement-limited. In units of the first
  null $\beta_0$ of the bistatic pattern (4 ... 8.5 deg here, about $\lambda/D$), the co-pol ratio
  and phase, a mixed co/cross response, entropy and $\bar\alpha$ hold to $0.6\ldots1.0\,\beta_0$
  (0.5 dB, 5°, 3 dB, 0.1, 10°), beyond the span's −3 dB point at $0.41\ldots0.50\,\beta_0$, and break
  down together at the null. The sphere does not depart anywhere up to 16°.
* Weak or monostatically zero channels move first: the trihedral's cross-pol when the
  bistatic plane is not its mirror plane (co/cross 3 dB down at $0.06\,\beta_0$, 20 dB at
  $0.3\,\beta_0$), and the residual co-pol of a 45° dihedral ($0.2\,\beta_0$). LDR-type descriptors
  are the sensitive ones.
* Symmetry decides whether and how fast. If the bistatic plane is a mirror plane of the target
  (sphere, trihedral roll 90°, dihedral 0°/90°), cross-pol never appears. A target with C2
  symmetry about the line of sight (any dihedral or plate) has $S(-\beta) = S(\beta)$, so its
  polarimetric change is second order (field $\propto\beta^2$: 12 dB per doubling of β). The
  trihedral (C3, no C2) has an odd cross-pol, first order (field $\propto\beta$: −49 / −39 / −34 dB
  re span at 0.25 / 1 / 2°, about 5 dB per doubling against the 6 dB of a pure first-order
  term). Checked with mirrored receivers (dihedral bit-identical; trihedral
  cross-pol sign-flipped to within its floor).
* Genuine and numerical cross-pol differ in convergence, not in smoothness. A genuine bistatic
  cross-pol keeps its level and phase as the ray density rises (trihedral at 1° and 2°: within
  0.7 dB and about 8° at 4, 8 and 16 rays/λ); the floor does not (monostatic −61 / −74 /
  −103 dB, phase unrelated between densities). The genuine cross–co phase is deterministic: for
  the trihedral it is linear in frequency, with the delay of its aperture rim relative to its
  corner (two-way path 2a/√3: 46 / 92 / 184 mm for edges a of 40 / 80 / 160 mm), −124 … −134°
  at 79 GHz; for the rolled dihedrals flat at 0/180°. The SBR+ floor is fixed relative to the
  monostatic span (−60 dB trihedral, −85 ... −100 dB sphere and dihedral at 4 rays/λ), so
  cross-pol is only defined well above it.
* Cloude–Pottier's stability is partly blindness: the trihedral's co/cross had fallen 20 dB by
  2.5° while H and $\bar\alpha$ stayed put to 6.5°.
* Aiming the receivers at the target changes the span (beam gain) but no polarimetric quantity:
  with ideal beams there is no geometric squint in the principal plane.
* Size and range (trihedral 40–160 mm at 1–16 m, dihedral 160 mm at 4 m): $\beta_0 \propto
  \lambda/D$ within 2 %, independent of range. The co-pol rule above holds in $\beta/\beta_0$ at
  every size and range, so the co-pol signature survives an array baseline $d$ while
  $d \lesssim 0.5\,R\lambda/D$ (at 79 GHz and 10 m with $d$ = 50 mm: targets below about 0.4 m).
* A trihedral's bistatic cross-pol depends on $\beta$, not on range, and falls about 9 dB per
  doubling of size: at 1° −39 / −43 / −57 dB re span for 40 / 80 / 160 mm. A calibration-size
  trihedral (≥ 10 cm) at array-scale angles stays near −50 dB, below real antenna cross-pol.
  Edge diffraction makes up most of it (physical optics alone: 10 dB lower), so the level needs
  a full-wave check.
* Large targets in their near field (160 mm dihedral at 4 m, far field at 27 m) have no lobe:
  the power ripples by ±3 dB across the retro-reflected footprint (+5 dB at 1°) while the
  co-pol ratio stays within 0.4 dB and the phase within 4°.

### Findings: supervisor's suggestions, more reflectors, split pairs (2026-10-21)

Same scene and conventions (now radar ones: LDR = cross-polar power over the co-polar power of
the same Tx polarization, Z_DR = |S_HH|²/|S_VV|²); "co/cross 3 dB down" above reads as LDR
3 dB up. Added: a slant trihedral (roll 45°), a 40 mm flat plate, and split antenna pairs.

* LDR is the most sensitive descriptor, but only where symmetry allows cross-pol: the
  trihedral at roll 0° and 45° (LDR −37 / −39 dB at 1°, rising 6 dB per doubling of β), never the
  sphere, the plate, the trihedral at roll 90° or the dihedral at 0°/90°. With real antennas
  (pattern cross-pol typically −20 … −30 dB) these bistatic levels would sit under the
  antenna's own LDR floor, except for small calibration-size trihedrals.
* The cross–co phase ∠(S_VH S_HH*) is deterministic wherever the cross-pol is genuine (trihedral:
  linear over the band with its rim-to-corner delay, −124 … −134° at 79 GHz); at the numerical
  floor it changes with the ray density and means nothing.
* Symmetry sets the order of every change (receivers displaced along H). A mirror plane x–z
  (perpendicular to the displacement) or C2 symmetry about the line of sight: co-pol changes to
  second order, cross-pol odd. A mirror plane in the bistatic plane: no cross-pol, but the
  co-pol may change linearly (trihedral roll 90°: Z_DR +0.2 dB, ∠(S_HH S_VV*) +1.1° at 1°).
  Neither (slant trihedral): both linear, so its Z_DR departs 0.5 dB at 3.5° against 7.5° at
  roll 0°.
* The flat plate behaves like the dihedral at roll 0° (same 40 mm aperture along H): first
  null near λ/L = 5.4°, Z_DR 0.5 dB at 4.75°, and no cross-pol at all (physical optics).
* Z_DR from split pairs (TxV –b– RxV –g– TxH –b– RxH, both co-pol pairs with the same pair
  angle β) agrees with a co-located dual-polarized pair to within 0.2 dB and 0.9° for β ≤ 2°
  with g = 0, for every target with a resolved co-pol (sphere, plate, dihedral 0°/22.5°,
  trihedral 0°/45°/90°; 2026-10-29: plate and slant trihedral added): at automotive geometry
  the split costs nothing with ideal beams. The dihedral at 45°, whose co-pol is 16 … 44 dB
  below its cross-pol and not reciprocal in SBR+, gives 0.9 dB and 9°: no Z_DR to speak of. It departs at larger β and gaps, most for the target without a mirror plane x–z or C2
  (trihedral roll 90°: −1.3 dB at 3° and −7.6 dB at 6° with g = 2b), and in spikes where the
  pairs' bisectors fall onto the target's monostatic aspect nulls (dihedral, g = 2b).
* Aiming the receivers changes only the span (a scalar gain fit leaves at most −48 dB inside
  the lobe), for ideal identical beams in a principal plane.

### Limitations (apply to every finding above)

* **Ideal antennas.** Parametric beams with identical H and V patterns, no pattern cross-pol,
  no phase-centre wander, no coupling, one phase centre per site. Every "aiming does not matter"
  and "the split does not matter" statement holds only for them: real patterns add the receive
  polarization matrix P_Rx(β) (E/H-plane gain ratio, pattern cross-pol, basis rotation off the
  principal planes), which affects monostatic off-boresight targets as much as bistatic ones.
* **Geometry.** Horizontal displacement only (bistatic plane x–y), targets on boresight at
  4 m (1–16 m for the trihedral scaling), canonical PEC targets 40 mm across (up to 160 mm),
  free space: no ground, multipath or clutter. Vertical displacement swaps the targets'
  symmetry roles and is not run.
* **Solver.** SBR+ is asymptotic (PO with PTD/UTD edges). Edge-dominated quantities, above all
  the trihedral's bistatic cross-pol, are model-dependent by about 10 dB (PO against PTD + UTD)
  and need a full-wave check. Cross-pol is defined only well above the numerical floor (about
  −60 dB re span for the trihedral, −85 … −100 dB otherwise, at 4 rays/λ).
* **Target near field.** The lobe rules need the target's far field (R > 2D²/λ); the 160 mm
  dihedral at 4 m has no lobe at all.
* **Processing.** Static frequency sweeps (74–84 GHz), frequency-domain beamforming towards a
  known target, no range gating, no Doppler, no calibration errors.

### Decisions (with their reasons)

| Date | Decision | Reason |
| --- | --- | --- |
| 2026-10-20 | Parametric beams, not the full-wave array, for the bistatic scenes | Isolate the target's bistatic physics from antenna effects; real patterns are a separate, one-way question |
| 2026-10-20 | Receivers on a ladder y_k = R tan(kβ), pooled by β from the exported positions | One solve gives three bistatic angles; β is never taken from a formula |
| 2026-10-20 | 4 rays per wavelength, sphere at mesh slider 9 | T4: converged above −55 dB re co-pol; T6: faceting at slider 5 |
| 2026-10-20 | Sphere, trihedral (0°, 90°), dihedral (0 … 90° in 22.5° steps) | User validation; 22.5° added for a mixed co/cross response |
| 2026-10-20 | Departure levels relative to each target's monostatic value, floors masked at 50 dB below the span, levels only in CSVs | Absolute levels assumed an odd-bounce target; figures stay unbiased |
| 2026-10-20 | Horizontal displacement only | Vertical noted as an extension |
| 2026-10-21 | Receiver aiming closed by the scalar fit (`reference/rx_aim/`) | User: sufficient proof; scope limited to ideal beams |
| 2026-10-21 | Scene frame ISO 8855 (x boresight) by default, `--boresight z` kept; exports record `target_steering_deg`, the analysis steers towards it | Consistency with the signal processing; the z scene's angles are spherical |
| 2026-10-21 | Radar conventions (Z_DR, LDR, ∠(S_HH S_VV*), ∠(S_cross S_co*)) under new keys | Supervisor's suggestion; new keys so that no old name changes meaning |
| 2026-10-21 | Figures split by signature group, labels in mathtext with power ratios as ratios on dB axes; summaries one figure per target family | User: the big figure was too big and " | S_HH | − | S_VV | " unclear; the old summary dropped series |
| 2026-10-21 | Slant trihedral (roll 45°) and flat plate added; slant dihedrals already covered | Supervisor: more reflectors |
| 2026-10-21 | Split pairs as mirror images about x–z, gap 0 and 2b, against a co-located pair of the same base | Supervisor's Z_DR layout; equal pair angles and gains isolate the split itself |
| 2026-10-24 | β (and β/β₀) as the swept and plotted variable, the design rule stated in d | Far-field scattering depends on d and R only through β; it fails in the target's near field (glossary: "Why β, not the spacing d") |

### Log entry

* 2026-10-20: Multi-receiver sweep done (sandbox branch `sandbox/bistaticity`, plan and log in
  `src/pomaar/simulator/investigation/investigation_log.md`): sphere, trihedral (2 rolls),
  dihedral (5 rolls) at 4 m, bistatic angles 0 ... 16°; six trust checks plus a β-sign symmetry
  check pass. Departure table and summaries in `data/processed/bistaticity/parametric_beams/summary/`;
  size and range scaling (trihedral 40–160 mm, 1–16 m) in `…/scaling/`; edge-diffraction check
  (PO, PTD, UTD).
  Extension noted: vertical receiver displacement.
* 2026-10-21: Pipeline validated end to end on a copy of the project (16 cases, 6 bugs fixed;
  `src/pomaar/simulator/investigation/validation/README.md`). Signatures switched to radar
  conventions (Z_DR, LDR, cross–co phase), figures split by group with clearer labels. Supervisor's
  suggestions run: slant trihedral, flat plate, and Z_DR from split antenna pairs
  (`…/pairs/pairs_zdr.png`): the split agrees with a co-located pair to 0.2 dB for β ≤ 2° with
  ideal beams. Limitations and a glossary of the aperture-physics terms in the log file.
* 2026-10-24: Presentation figures in `reports/bistaticity-analysis/` (parametric beams and the
  horn cluster re-analysed). Corrected the trihedral's cross–co phase: linear over the band with
  its rim-to-corner delay 2a/√3, not flat at +135°; a floor's phase is smooth too and is told
  apart by its change with ray density.

## Glossary: geometry, aperture physics and descriptors

### Geometry

* **Bistatic angle β.** The angle at the target between the directions to the Tx and to the Rx
  (`bistatic_angle_deg`); monostatic: β = 0. For a Tx and Rx at ±d/2 about a centre and a
  target at range R on the perpendicular bisector, tan(β/2) = d/(2R). In the multi-receiver
  ladder the range is measured from the Tx, and Rx k sits at y = R tan(kβ), so tan(kβ) = y/R
  exactly; both are β ≈ d/R to first order. The analysis always takes β from the exported
  positions, never from a formula.
* **Why β, not the spacing d, is the swept variable.** In the target's far field
  (R > 2D²/λ), the incident and scattered waves are plane over the target, so S depends only
  on the incident and scattering directions, and d and R enter only through β ≈ d/R. One β
  sweep therefore covers every (d, R): the scaling run gave the same β₀ and the same trihedral
  cross-pol at 1, 4 and 16 m for equal β. A d sweep at one range would hide the R and D in the
  rule. The design rule follows from β ≲ 0.5 β₀ and β₀ ≈ λ/D as d ≲ 0.5 Rλ/D. The binding
  case is the shortest range and the largest target, d_max = 0.5 R_min λ/D_max. At the
  far-field boundary R = 2D²/λ it gives d_max ≈ D.
  * Near-field breakdown: for R < 2D²/λ the wave front curves across the target, the
    retro-reflected beam has a geometric footprint with Fresnel ripple instead of a lobe
    (160 mm dihedral at 4 m: ±3 dB across 0–2.5°), and d and R act separately, so β alone no
    longer describes the scene.
* **Bistatic plane.** The plane through Tx, Rx and target: x–y for the horizontal displacement
  studied here, x–z for a vertical one.
* **Bistatic bisector.** The direction halfway between the Tx and Rx directions, seen from the
  target. For small β and smooth, perfectly conducting targets, a bistatic measurement behaves
  like a monostatic one along the bisector at the frequency scaled by cos(β/2) (the
  monostatic–bistatic equivalence theorem). Multiple bounces, shadowing and edges break it,
  which is why the trihedral and dihedral depart from their monostatic signatures at all.
* **Split pairs.** Two separate bistatic pairs, one per co-pol channel (TxV–RxV and TxH–RxH,
  base b each, gap g between them). Their bisectors lie (b + g)/R apart, so the two channels see
  the target from slightly different aspects.
* **Mirror plane and C2 symmetry (of a target).** A mirror plane is a plane whose reflection
  leaves the target unchanged; C2 is a 180° rotation about the line of sight that does so. Both
  belong to the target, not to the geometry. With M = diag(−1, 1) (H flips, V stays) and the
  receivers displaced along H:
  * a mirror plane x–z gives S(−β) = M S(β) M: co-pol even in β (second order), cross-pol odd;
  * a mirror plane in the bistatic plane (x–y) forces the cross-pol to zero at every β, but
    leaves the co-pol free to change linearly;
  * C2 gives S(−β) = S(β): everything even, so all changes are second order.
* **Principal plane.** A plane containing the boresight and the H or V axis (x–y or x–z here).
  In it, the antennas' H/V basis is not rotated relative to the scattering basis.
* **Polarization basis.** H along y and V along z (ISO 8855 scene frame). S[rx, tx] with 0 = H
  and 1 = V (back-scatter alignment), so S_HV is H-receive, V-transmit.

### Aperture physics

* **Bistatic pattern.** A target's scattering as a function of the receive direction for a
  fixed transmit direction: the bistatic RCS σ(β). It is a property of the target, not of the
  antennas.
* **Why the plotted span is the bistatic pattern.** The span ‖S‖² = Σ|S_pq|² is the
  polarization-summed received power. With the receivers aimed at the target (rxAim 1), the
  antenna gains are the same at every β, and the paths change by only R/cos(kβ) (0.34 dB at
  16°). So the span against β is the bistatic pattern up to a constant (the radar equation's
  gains and ranges), summed over the polarizations. With the receivers on boresight (rxAim 0),
  the receive beam's roll-off is added on top: 3(β/30°)² dB for the 60° beams.
* **Why it has lobes.** A flat or corner reflector re-radiates like an aperture antenna. In
  physical optics (PO), the incident field induces surface currents J = 2 n̂ × H_inc on the lit
  faces. Their phases make them add up in the retro (or specular) direction. Moving the receiver
  by β adds a path difference that grows linearly across the aperture, so
  σ(β) ∝ |∫_A e^{jk(r̂_inc + r̂_s)·r′} dA|², a Fourier transform of the effective aperture. That
  gives a sinc-like **main (retro) lobe** around the monostatic direction, then **nulls** and
  **sidelobes**.
* **First null β₀ ≈ λ/D**, D the effective aperture along the bistatic plane:
  * trihedral, 40 mm edges: 7.8° (D ≈ 28 mm);
  * plate, and dihedral at roll 0 (40 mm along H): λ/L = 5.4°, seen at 5.75°;
  * dihedral at roll 90 (retro-reflecting across the seam, D = 2L cos 45°): 3.85°, seen at
    4.0°.

  β₀ scales as 1/D and does not depend on range (checked over 40–160 mm and 1–16 m). The
  span's −3 dB point lies at 0.41–0.50 β₀.
* **Monostatic aspect pattern.** The monostatic RCS against the target's aspect angle θ. Its
  phase gradient is twice the bistatic one (2k sin θ across the aperture), so its first null
  lies at λ/(2D), half of β₀ (2.7° for 40 mm). It matters for the split pairs, whose bisectors
  sit off boresight.
* **Lobe-limited.** Descriptors of strong channels (Z_DR, co-pol phase, H, ᾱ) hold to
  0.6–1.0 β₀ and change together near the first null, where the dominant mechanism vanishes and
  weaker ones (edges, single or double bounces) compete. In displacement terms, d ≲ 0.5 Rλ/D.
* **Far field of the target.** The lobe structure needs R > 2D²/λ (1.7 m for a 40 mm trihedral,
  27 m for a 160 mm dihedral at 79 GHz). Closer in, the retro-reflected beam has a geometric
  footprint with Fresnel ripple instead of a lobe.
* **Sphere.** Smooth and doubly curved: σ ≈ πa² for small β and no lobe. The reference for a
  bistatically invariant target.
* **Edge diffraction (PTD/UTD).** PO alone ignores the edges. PTD (physical theory of
  diffraction) adds edge currents, and UTD (uniform theory of diffraction) adds diffracted
  rays. About 10 dB of the trihedral's bistatic cross-pol comes from edges (PO −49 dB against
  PTD + UTD −39 dB re span at 1°). Its −64 dB monostatic "floor" is UTD-ray noise.
* **Numerical floor.** The level below which SBR+ output is ray-sampling noise. It is fixed
  relative to the monostatic span (−60 dB trihedral, −85 … −100 dB sphere and dihedral, at
  4 rays/λ), and it falls with ray density. Below it, ratios are meaningless and phases random.

### Descriptors (power ratios in dB, phases in deg)

* **Span** ‖S‖² = Σ|S_pq|²: with aimed receivers, the bistatic pattern.
* **Differential reflectivity** Z_DR = |S_HH|²/|S_VV|² and the **co-pol phase difference**
  ∠(S_HH S_VV*):
  * 0 dB and 0° for odd bounces (sphere, trihedral, plate);
  * 0 dB and 180° for a dihedral at roll 0.
* **Linear depolarization ratio** LDR_H = |S_VH|²/|S_HH|² (Tx H) and LDR_V = |S_HV|²/|S_VV|²:
  the cross-polar power normalized to the co-polar power of the same transmit polarization.
* **Cross–co phase** ∠(S_VH S_HH*) and ∠(S_HV S_VV*): meaningful only where the cross-pol stands
  above its floor.
* **Cross-pol ratio** |S_HV|²/|S_VH|² and the **cross-pol phase difference** ∠(S_HV S_VH*).
  The **normalized asymmetry** is |S_HV − S_VH|²/‖S‖². Monostatically, reciprocity fixes these
  to 0 dB, 0° and zero; bistatically they measure how far S departs from symmetry.
* **Bistatic Pauli fractions** |k_i|²/‖k‖², with k = (HH+VV, HH−VV, HV+VH, j(HV−VH))/√2:
  * k1: odd bounce;
  * k2: even bounce;
  * k3: 45°-rotated even bounce;
  * k4: the antisymmetric part, which only a bistatic geometry (or a non-reciprocal medium)
    produces.
* **Cloude–Pottier H, A, ᾱ** of the 4×4 band-averaged coherency matrix:
  * entropy H: 0 = one mechanism;
  * anisotropy A;
  * mean alpha ᾱ: 0° odd bounce, 45° dipole, 90° even bounce.

  They summarize, and so they can be blind to changes in weak channels: the trihedral's LDR rose
  20 dB while H and ᾱ stayed put.
* **Departure angle.** The first β at which a descriptor moves from its monostatic value by a
  stated level (`summary_departures_aim<a>.csv`). It is an interpretation, so it appears only in
  CSVs, never in figures.
* **Scale invariance.** Every descriptor except the span is a ratio, a phase difference or a
  span-normalized quantity, so S → gS (one complex scalar) leaves it unchanged. That is why
  receiver aiming, a common gain with ideal beams, changes only the span.
* **Single-pol / dual-pol (array topology).** A single-pol array uses separate
  single-polarized elements per polarization state, so the polarization channels have different
  phase centres; a dual-pol array uses dual-polarized elements, so they share one (`spacingHV` =
  0 here). The split pairs are the bistatic analogue of single-pol. Formerly type 1 / type 2;
  with their processing they are the research proposal's architecture A (single-pol,
  power-only) and B (dual-pol, coherent).
