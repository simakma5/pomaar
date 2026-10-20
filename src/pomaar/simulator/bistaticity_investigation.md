# Bistaticity investigation: running plan and log

Working notes of the sandboxed investigation behind the task note "Bistaticity analysis"
(Obsidian). A distillation goes to that note's Findings and Log; this file keeps the details.

## Sandbox

| What | Where |
| --- | --- |
| Branch | `sandbox/bistaticity` (git worktree `~/Repositories/pomaar-sandbox`; `main` and its working tree untouched) |
| Baseline | `9dfe44b`: the uncommitted simulator state of `main` on 2026-10-20 |
| AEDT project | `~/Repositories/pomaar/data/interim/sandbox/bistaticity/bistaticity_sandbox.aedt`, a copy of `~/Projects/AEDT/bistaticity-analysis/bistaticity_analysis.aedt` (12:21), open in the GUI session next to the original |
| Scripts | `src/pomaar/simulator/investigation/` (run with `ansys-vnc exec python ...` from the worktree, so its `.venv` and code are used); they refuse any project but the sandbox copy |
| Data | `~/Repositories/pomaar/data/interim/sandbox/bistaticity/` (exports, figures) |

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
* `spacingHV` (H at z = +s/2, V at z = -s/2 in every site) stays 0 here (type 2): the
  polarization channels share each site's phase centre.
* Names: instances `Tx_H`, `Tx_V`, `Rx_<k>_H`, `Rx_<k>_V`; SBR+ appends `_p1` to their ports,
  which `parse_port` now strips.

### Plan

- [x] M0 Sandbox: worktree and branch, project copy, plan and log.
- [x] M1 Code: `parse_port` accepts SBR+ native-antenna ports (`<component><index>_p<N>`);
  sandbox helpers (`investigation/`).
- [x] M2 Build the multi-receiver scene in the sandbox project (`build_scene.py`;
  `MultiRxTrihedral`: 40 mm trihedral at 4 m, bistaticAngle 1 deg, rxAim 1, spacingHV 0).
- [x] M3 Trust checks (small solves), then **stop for validation** (all pass; see the log):
  - T1 level: monostatic trihedral |S_HH| = |S_VV| against the radar equation;
  - T2 polarization labels: relative to a trihedral, a dihedral with its seam along H flips
    S_VV and keeps S_HH (the planned PEC strip along H turned out uninformative);
  - T3 phase: arg S follows -k (r_Tx + r_Rx) of each antenna's own position;
  - T4 numerical cross-pol floor against ray density (4, 8, 16 rays per wavelength);
  - T5 the co-located monostatic receiver (Rx_0 on top of Tx) solves cleanly;
  - T4b symmetry: the trihedral's cross-pol in the bistatic plane vanishes when that plane is
    its mirror plane (targetRoll 90 deg) and may not otherwise (targetRoll 0).
- [x] M4 Production sweep (validated 2026-10-20, with a dihedral at 22.5 deg added):
  - designs (all 40 mm across, at 4 m): `MultiRxSphere` (radius 20 mm, mesh slider 9, T6),
    `MultiRxTrihedral` at targetRoll 0 and 90 deg, `MultiRxDihedral` at targetRoll 0, 22.5,
    45, 67.5 and 90 deg (a linear roll sweep; 67.5 mirrors 22.5);
  - `bistaticAngle` 0.25 ... 8 deg in 0.25 deg steps (32 values; with Rx_2 the bistatic
    angles reach 16 deg), each with rxAim 0 and 1: 512 solves, about 3.5 h
    (`investigation/run_production.sh`, log `production.log` in the sandbox data folder);
  - exports: one per (design, rxAim, targetRoll) over `bistaticAngle`, with the scene
    geometry (`positions_m`, `target_m`) in the metadata (`export_sweep.py`);
  - analysis: `bistaticity_analysis --per-receiver` pools every (receiver, variation) by its
    bistatic angle (`pool_receivers`); a summary script compares targets and rxAim.
  - [x] T6 sphere faceting: at the default mesh (slider 5) the 20 mm sphere was 2.4 dB above
    the radar equation with +-2.5 dB ripple over the band, and its co-pol ratio drifted
    0.18 dB / 1.7 deg at 4 deg bistatic; at slider 8 and 9: -113.1 dB against -113.3 dB,
    +-0.6 dB ripple, 0.027 dB / 0.43 deg at 4 deg (8 and 9 agree), cross-pol floor -90 / -95
    dB. The sphere runs at slider 9; the plate targets are not faceted by curvature.
  - [x] T7 symmetry in the sign of the bistatic angle (receivers mirrored to -y, 3 solves): a
    target with C2 symmetry about the line of sight gives S(-beta) = S(beta), one with a
    mirror plane x-z gives S(-beta) = M S(beta) M with M = diag(-1, 1) (H flips).
- [x] M5 Summary of the action point for the task note (Findings, Log): see "Distillation".

### Extension action points (not in this round)

- [ ] Vertical receiver displacement (bistatic plane x-z, along V), next to the horizontal
  one (x-y, along H) studied here: the trihedral's and dihedral's symmetry planes swap roles.

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
  sits at the floor. The span falls 0.4 dB by 16 deg with aimed receivers (the longer receive
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
  band and over 0.25 ... 6 deg) and HV = VH within 1 dB up to about 4 deg (reciprocity holds; the
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
  'undefined' (and is left out of `summary_change_*.png` per design).
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
  cross-pol at roll 0 grows 17 dB per decade (-49 / -44 / -39 / -34 dB re the monostatic span
  at 0.25 / 0.5 / 1 / 2 deg): field proportional to beta. The rolled dihedrals' emerging
  channels grow 40 ... 47 dB per decade: k2 at 45 deg (-82 / -68 / -55 dB at 0.5 / 1 / 2 deg)
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
  deg, flat over the band; rolled dihedrals: 0 / 180 deg). Cloude-Pottier's stability holds
  (H and alpha move only at 0.8 ... 1.0 beta_0) but is partly blindness: at roll 0 the
  trihedral's co/cross had fallen 20 dB by 2.5 deg while H and alpha stayed put to 6.5 deg.

## Distillation (for the task note)

### Findings: multiple receivers, boresight canonical targets at 4 m (2026-10-20)

Scene: ideal parametric beams (60 deg, type 2: H and V share each site), one Tx and three Rx
sites at bistatic angles 0, beta, 2 beta (beta = 0.25 ... 8 deg, so 0 ... 16 deg in all),
horizontal displacement only; sphere, trihedral (roll 0, 90 deg) and dihedral (roll 0, 22.5,
45, 67.5, 90 deg), all 40 mm across. 512 SBR+ solves, sandbox branch `sandbox/bistaticity`.

- Strong-channel descriptors are lobe-limited, not displacement-limited. In units of the first
  null $\beta_0$ of the bistatic pattern (4 ... 8.5 deg here, about $\lambda/D$), the co-pol ratio
  and phase, a mixed co/cross response, entropy and $\bar\alpha$ hold to $0.6\ldots1.0\,\beta_0$
  (0.5 dB, 5°, 3 dB, 0.1, 10°), beyond the span's −3 dB point at $0.41\ldots0.50\,\beta_0$, and break
  down together at the null. The sphere does not depart anywhere up to 16°.
- Weak or monostatically zero channels move first: the trihedral's cross-pol when the
  bistatic plane is not its mirror plane (co/cross 3 dB down at $0.06\,\beta_0$, 20 dB at
  $0.3\,\beta_0$), and the residual co-pol of a 45° dihedral ($0.2\,\beta_0$). LDR-type descriptors
  are the sensitive ones.
- Symmetry decides whether and how fast. If the bistatic plane is a mirror plane of the target
  (sphere, trihedral roll 90°, dihedral 0°/90°), cross-pol never appears. A target with C2
  symmetry about the line of sight (any dihedral or plate) has $S(-\beta) = S(\beta)$, so its
  polarimetric change is second order (field $\propto\beta^2$, 40 dB/decade). The trihedral (C3, no
  C2) has an odd cross-pol, first order (field $\propto\beta$, 20 dB/decade: −49 / −39 / −34 dB re
  span at 0.25 / 1 / 2°). Checked with mirrored receivers (dihedral bit-identical; trihedral
  cross-pol sign-flipped to within its floor).
- The erratic cross-pol phase is the numerical floor: genuine bistatic cross-pol has a
  deterministic phase (trihedral about +135°, flat over the band; rolled dihedrals 0/180°).
  The SBR+ floor is fixed relative to the monostatic span (−60 dB trihedral, −85 ... −100 dB
  sphere and dihedral at 4 rays/λ), so cross-pol is only defined well above it.
- Cloude–Pottier's stability is partly blindness: the trihedral's co/cross had fallen 20 dB by
  2.5° while H and $\bar\alpha$ stayed put to 6.5°.
- Aiming the receivers at the target changes the span (beam gain) but no polarimetric quantity:
  with ideal beams there is no geometric squint in the principal plane.
- Implication for type 1 vs type 2 (hypothesis, to be tested with target size and range): an
  array-scale displacement $d$ of a few λ at metres puts $\beta/\beta_0 \approx dD/(R\lambda)$ near
  $10^{-2}$, where strong channels are untouched and a trihedral's bistatic cross-pol is about
  −60 dB re span. It grows 20 dB per decade of $dD/(R\lambda)$, so it matters for large or
  close targets and for weak-channel descriptors first.

### Log entry

- 2026-10-20: Multi-receiver sweep done (sandbox branch `sandbox/bistaticity`, plan and log in
  `src/pomaar/simulator/bistaticity_investigation.md`): sphere, trihedral (2 rolls), dihedral
  (5 rolls) at 4 m, bistatic angles 0 ... 16°; six trust checks plus a β-sign symmetry check
  pass. Departure table and summaries in `data/interim/sandbox/bistaticity/summary/`.
  Extension noted: vertical receiver displacement.
