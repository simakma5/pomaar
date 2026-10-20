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
  2 deg) must agree: a built-in consistency check.
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
- [ ] M4 Production sweep (after validation), proposed:
  - designs: sphere (radius 20 mm), trihedral at targetRoll 0 and 90 deg (bistatic plane off
    and on its mirror plane), dihedral at targetRoll 0, 45 and 90 deg; all 40 mm at 4 m;
  - `bistaticAngle` 0.25 ... 8 deg in 0.25 deg steps (32 values; with Rx_2 the bistatic
    angles reach 16 deg), each with rxAim 1 and 0: 64 solves of about 25 s per design, under
    half an hour each;
  - analysis: a per-receiver mode for `bistaticity_analysis` that pools every (variation,
    receiver) into signatures against the bistatic angle k bistaticAngle, the co/cross-pol
    signatures (Z_DR, LDR analogues) included; equal angles from different receivers as the
    consistency check; rxAim 1 against 0 for the geometric squint.
- [ ] M5 Summary of the action point into the task note (Findings, Log).

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
