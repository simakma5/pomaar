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
- [ ] M1 Code: `parse_port` accepts SBR+ native-antenna ports (`_p<N>`); sandbox helpers.
- [ ] M2 Build the multi-receiver scene in the sandbox project (script, reproducible).
- [ ] M3 Trust checks (small solves), then **stop for validation**:
  - T1 level: monostatic trihedral |S_HH| = |S_VV| against the radar equation;
  - T2 polarization labels: a thin PEC strip along y (H) gives |S_HH| >> |S_VV|;
  - T3 phase: arg S follows -k (r_Tx + r_Rx) of each antenna's own position;
  - T4 numerical cross-pol floor against ray density (4, 8, 16 rays per wavelength);
  - T5 the co-located monostatic receiver (Rx_0 on top of Tx) solves cleanly.
- [ ] M4 Production sweep of `bistaticAngle` (after validation) for sphere, trihedral and
  dihedral (roll 0, 45, 90 deg), then the analysis: signatures against the bistatic angle.

## Log

### 2026-10-20

* M0: created the worktree and branch `sandbox/bistaticity` from `main` (2afa6be) and
  committed the uncommitted simulator state of `main` as the baseline (9dfe44b). Copied the
  project (saved 12:21) to the sandbox as `bistaticity_sandbox.aedt`: AEDT cannot open two
  projects of the same name, and the distinct name keeps every script away from the original.
