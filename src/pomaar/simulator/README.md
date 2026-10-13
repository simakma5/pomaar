# Polarimetric MIMO simulation sub-project

This sub-project automates a three-stage simulation pipeline linking **Ansys HFSS (full-wave)**, **Ansys HFSS SBR+ (ray-tracing)**, and **Python DSP post-processing** for polarimetric MIMO automotive radar array analysis.

---

## 1. Prerequisites and environment setup

The pipeline requires **Ansys Electronics Desktop (AEDT) 2025.2 or newer** and a Python environment with PyAEDT (`pyaedt`), provided by this repository's `uv` environment (`uv sync`).

### Option A: Running in the AEDT container (Linux, recommended)

On the Fedora workstation, AEDT runs in the `ansys-vnc` container (see `~/Repositories/ansys-vnc`). Python scripts must run **inside** that container: from the host, PyAEDT cannot see the container's gRPC session and would silently start a native AEDT instead (the builder refuses to do so). The container mounts `~/Repositories` and the host's `uv` Python installations, so the project's `.venv` works there unchanged:

```bash
ansys-vnc start                      # once; then open the GUI from the VNC desktop or `ansys-vnc aedt`
cd ~/Repositories/pomaar
ansys-vnc exec hfss_array_builder <path_to_project.aedt> <unit_cell_design_name> [layout.yaml]
```

To run the scripts from any directory, link them once: `ansys-vnc link --host bistaticity_analysis ~/Repositories/pomaar` writes wrappers into `~/.local/bin`. `hfss_array_builder`, `sbr_simulator` and `sbr_results` then run in the container, and `bistaticity_analysis` runs on the host. Rerun the command after adding console scripts.
`ansys-vnc exec` keeps the current directory, pins the process to the compute cores and puts the nearest `.venv` (here the project's) first on `PATH`, so the console scripts and `python` run by name. The builder connects to the running GUI session's gRPC port (found automatically; `--port` overrides it) and starts a non-graphical AEDT in the container if none is listening.

### Option B: Running natively (Windows / supported Linux)

Ensure that:

1. `ANSYSEM_ROOT252` (or your corresponding AEDT installation path) is set in your environment variables.
2. The `pomaar` package dependencies are installed (`uv sync`).
3. On Linux outside a container, set `POMAAR_ALLOW_NATIVE_AEDT=1` to let the builder start or connect to AEDT.

---

## 2. Running the pipeline

### Module 1: HFSS full-wave array synthesis (`hfss_array_builder.py`)

This script automates the creation of a planar MIMO array on a single contiguous PCB board by copying, replicating, and boolean-cutting template geometries from an isolated unit-cell element design.

To run (prefix with `ansys-vnc exec` in the container setup):

```bash
hfss_array_builder <path_to_project.aedt> <unit_cell_design_name> [layout.yaml] [--unit-cell POL=DESIGN ...]
```

* **Example:**

  ```bash
  hfss_array_builder ~/Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt "ApertureCoupledPatch"
  ```

* **Polarization:** an element's polarization comes from its unit-cell template, not from its placement. Port sheets carry their native polarization in the name (`PortSheet_H`, `PortSheet_V`); a unit cell has one such port (single-pol) or two (dual-port). Further unit cells are given per polarization set, `--unit-cell h=PatchH`, `--unit-cell v=PatchV` for single-pol elements and `--unit-cell hv=DualPatch` for dual-polarized ones (`polarization: hv`), repeatable, or `metadata.unit_cells: {h: PatchH, hv: DualPatch}` in the layout, and a layout element may name its own with `unit_cell:`. Elements without one use the source design, turned by their `yaw`. Either way the builder checks that R(yaw) applied to each port's native E-field lies along x for `polarization: h` and y for `v` (`hv` needs a dual-port cell) and aborts otherwise, before it touches any design. Untagged legacy sheets (`PortSheet`) take `--native-polarization h|v` (source design) or their `--unit-cell` key.
* **Unused ports:** an element that uses only one port of a dual-port cell keeps the other feed, terminated in a matched wave port `Port_<label>_<P>_Load`, so that it does not become an open stub; `--drop-unused-ports` removes it instead. `sbr_simulator` skips `_Load` ports.
* **Field signs:** turning an element by 180° flips its E-field. After building, the builder prints every port's polarization and sign (+1 along +x for H, +y for V) and writes them, with the phase centres, to `<TargetDesign>_ports.json` next to the project. The measured matrix is diag(s_RxH, s_RxV) S diag(s_TxH, s_TxV), so a role whose H and V ports have opposite signs (pair sign −1) negates the cross-pol channels against the co-pol ones; a sphere calibration cannot reveal this, and the builder warns about it.
* **What it does:** Creates/updates a target design named `ApertureCoupledPatchMimoArray`, syncs design variables, sizes the PCB board to perfectly fit feedlines, duplicates Tx and Rx elements according to grid positions, pre-enables post-processing, assigns wave ports (lumped-port fallback), excites the first Tx port for far-field reports, draws a $\lambda_0/4$ vacuum radiation Airbox, and sets up Theta-Phi far-field spheres in `PhaseCentreCS` with realized-gain pattern and directivity-based XPD reports. The far-field phase reference is set per port to each element's CS (`--phase-reference-mode per-port`, the default), so that a linked SBR+ antenna launches every port's rays from its own element phase centre; `single` references all ports to the array centroid `PhaseCentreCS` instead, which is only valid for targets in the far field of the whole array.

#### Layouts (`layouts/`)

Each YAML is regenerated by the `generate_<topology>.py` next to it (`--help` lists the options). All
assume a V-polarized source design whose feed points to −y at yaw 0.

| Layout | Elements | Status |
| --- | --- | --- |
| `contiguous_ula_3x4`, `contiguous_ula_4x4` | V-only Rx and Tx ULAs | Usable; single-pol, so no polarimetric quad |
| `spatial_diversity_2d` | V-only L-shaped Rx and Tx row | Usable; single-pol, so no polarimetric quad |
| `polarimetric_cluster` | Tx/Rx × H/V, H by turning the V cell (feeds pointing outwards) | Usable; pair signs +1, one quad `Rx/Tx` |
| `polarimetric_ula` | Dual-polarized `hv` Rx elements, alternating V/H Tx | Needs a dual-port cell, `--unit-cell hv=<Design>` (or `generate_polarimetric_ula.py --dual-pol-cell <Design>`); two Tx per polarization, so the processor needs the four ports of a quad explicitly |

The earlier `polarimetric_ula` placed a V and an H copy of the single-pol cell at every Rx position. That geometry overlaps and gave both roles a pair sign of −1, so it is no longer generated.

### Module 2: SBR+ target solver (`sbr_simulator.py`)

This script builds an SBR+ design next to a solved HFSS antenna design of the same project: it places a canonical PEC target, links the HFSS design as the antenna source, pairs every Tx port with every Rx port, copies the frequency sweep, creates quick-look S-parameter reports and optionally solves.

To run (prefix with `ansys-vnc exec` in the container setup; `uv sync` installs the script):

```bash
sbr_simulator <path_to_project.aedt> <hfss_design_name> [--target dihedral] [--size 40mm] [--distance 4m] [--map-variables VAR ...] [--field-type farfield] [--simulate | --no-simulate]
```

* **Example:**

  ```bash
  sbr_simulator ~/Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt DualPolHornCluster --target sphere --map-variables copolarSpacingLambda
  ```

* **Units:** `--distance` is in metres and `--size` in mm unless a unit is given (`4m`, `4 metres`, `4 meters` and `400mm` all work); metres are written to AEDT as its unit token `meter`.
* **What it does**, in the design `Linked<hfss_design_name>` (`--design` overrides; an existing one is only replaced after confirmation or `--overwrite`):
  1. Creates `AntennaCS` at the global origin and places the linked antenna in it (boresight +z).
  2. Creates `TargetCS` in `AntennaCS` at the spherical position `X = targetDistance cos(targetAzimuth) sin(targetElevation)`, `Y = targetDistance sin(targetAzimuth) sin(targetElevation)`, `Z = targetDistance cos(targetElevation)`; `targetElevation` is therefore the polar angle off boresight. The CS's Euler angles ZYZ are `(Phi, Theta, Psi) = (targetAzimuth, targetElevation, targetRoll - targetAzimuth)`: the target's local z axis points radially away from the antenna without twisting, and `targetRoll` is a pure rotation about the line of sight (e.g. a 45° dihedral). Edit the TargetCS angles directly for other orientations.
  3. Creates the target in `TargetCS`, facing the antenna along local -z, sized by `targetSize` (`--size`) and made PEC by a Perfect E boundary (sheets carry no material); its SBR+ surface roughness and height standard deviation are 0 unless `--roughness` / `--height-deviation` are given:
     * sphere: PEC solid of radius `targetSize`, centred at `TargetCS`;
     * dihedral: two square `targetSize` × `targetSize` sheets, i.e. `targetSize` is the seam length and the sheet depth. The seam lies on the local x axis through `TargetCS` (parallel to H at zero roll), which is also the phase reference of the double bounce. One sheet is drawn, tilted by 45° and copied by a 180° rotation about the symmetry axis (local z);
     * trihedral: three right-isosceles triangular sheets with edge length `targetSize`, the corner at `TargetCS` and the symmetry axis pointing at the antenna. One sheet is drawn and copied by 120° rotations about local z.

     The defaults (40 mm at 4 m) keep all three in their own far field at 79 GHz (2D²/λ ≤ 3.4 m); the script warns when the chosen range is shorter. `--target none` builds the same scene without a target; its solve is the background (the direct Tx→Rx coupling between the proxy ports), a cross-check of the incident solution below (use `--design` to keep both designs).
  4. Links the HFSS design (solution `--source-solution`, default its first frequency sweep) with *Simulate source design as needed*, *Model Visualization* of all its objects except the radiation boundary objects (`RadiatingSurface`, the builder's `Airbox`), and far-field (standard) or near-field (bumper/radome integration) results, asked interactively unless `--field-type` is given. `--map-variables` maps independent HFSS design variables to SBR+ design variables of the same name, which are created from the HFSS value when missing, so they can be swept in the SBR+ design.
  5. Selects Tx/Rx: ports whose names contain `Tx`/`Rx` and end in `H`/`V` (e.g. `TxH`, `Port_Tx_1_V`); every Tx transmits to every Rx. The builder's matched `_Load` ports are skipped.
  6. Creates `Setup` (4 rays per wavelength, 5 bounces, PTD + UTD, creeping waves) whose sweep `Sweep` copies the linked HFSS sweep.
  7. Creates quick-look reports. SBR+ splits the solution `Setup : Sweep` (total) into `Setup : Sweep_Incident`, the direct Tx→Rx coupling between the proxy ports, and `Setup : Sweep_Scattered`, the target return (total = incident + scattered). The coupling dominates the total at short Tx/Rx spacings, so target signatures must be read from `Setup : Sweep_Scattered`. Reports: `S-parameters dB` (total) and `S-parameters dB (scattered)` for all Tx/Rx pairs, and, per polarimetric quad from the scattered solution, `Co-pol VV/HH` (`S_VV/S_HH` in dB and degrees; ideally 0 dB for both, with the dihedral's phase 180° away from the sphere's) and `Cross-pol isolation` (`S_VH/S_HH` for Tx H, `S_HV/S_VV` for Tx V, in dB). A quad pairs a dual-polarized Rx with a dual-polarized Tx: an element with an H and a V port of the same name prefix (`Port_Rx_1_H`/`_V`), or a role's single H and V ports when it has exactly one of each (the polarimetric cluster). When the HFSS design comes with the builder's port table (`<design>_ports.json` next to the project), the VV/HH phase is corrected for the port field signs.
  8. Asks whether to solve (`--simulate` / `--no-simulate` decide silently). Solving also solves the HFSS design where its linked solution is missing.

  `SbrSimulationManager.export_results(path, sweep_variable)` writes the solved scattered solution to an `.npz` for Module 3 (via `sbr_results.py`, which is also a standalone CLI; see §4.2).

### Module 3: Offline polarimetric post-processing (`polarimetry_processor.py`, `bistaticity_analysis.py`)

Pure NumPy, no AEDT: `MimoPolarimetryProcessor` loads an exported sweep and computes polarimetric signatures; `bistaticity_analysis` plots them against the Tx–Rx spacing. Section 4 documents the full chain, its validation and the first results.

### Spelling

The code uses British spelling with Oxford -ize (centre, metre, polarization). Command-line flags also accept the American spelling (`--center-freq`, `--use-existing-phase-center`), layouts may use `centre_frequency_ghz` or the older `center_frequency_ghz`, and names read from AEDT designs are looked up in both spellings (e.g. a unit cell's `PhaseCentreCS` or `PhaseCenterCS`, `--map-variables centreWavelength` for an HFSS `centerWavelength`). Lengths read from AEDT accept `meter(s)` and `metre(s)`.

---

## 3. Unit-cell design rules (strict modelling constraints)

To ensure the automated builder script runs successfully without manual alignment, any new unit-cell HFSS design **must** follow these strict modelling conventions:

> [!IMPORTANT]
> **Unit-cell modelling rules:**
>
> 1. **PCB layer naming:**
>    * Substrates and ground planes must follow the pattern `L{Order}_{Type}`:
>      * Substrates: e.g. `L12_Substrate`, `L23_Substrate` (must end with `_Substrate`).
>      * Ground planes: e.g. `L2_Ground` (must end with `_Ground`).
> 2. **Active geometries and ports:**
>    * Copper patches/feedlines must have unique layer names (e.g. `L1_Patch`, `L3_Trace`).
>    * Port excitation faces must be named `PortSheet` or `PortSheet{N}`, followed by their native polarization: `PortSheet_V`, `PortSheet1_H`. A dual-port cell has one of each; route both feeds to board edges that stay on the outside of the array, since wave ports must lie on the radiation boundary.
> 3. **Boolean dummy solids:**
>    * Slots, holes, or feed cutouts must be modelled as vacuum solids and named `f"{operation}_{target}"`:
>      * Example: `Subtract_L2_Ground` will automatically be replicated at every grid position and subtracted from the ground plane.
> 4. **`PhaseCentreCS` coordinate system:**
>    * The unit cell should contain a relative coordinate system named exactly **`PhaseCentreCS`** at the radiation phase centre of the unit cell.
>    * The script automatically retrieves `PhaseCentreCS` offsets to rotate/shift elements during replication. If it is missing (or you decline to reuse it), the builder fits it by weighted least squares to the co-polar far-field phase of the solved unit cell (solving the unit-cell setup first if needed); declining the fit falls back to a `[0,0,0]` offset.
>    * A dual-port cell may add `PhaseCentreCS_H` / `PhaseCentreCS_V` for its two ports (each port's far-field phase reference); it is never fitted.
> 5. **Several unit cells in one array:** all of them must have the same board stack (layer names, materials and z extents) and must not define a design variable differently (rename such variables per unit cell). The analysis setup comes from the source design, and the board stack is placed at the source design's phase-centre height.

---

## 4. Bistaticity study: polarimetric signatures from SBR+ sweeps

How much Tx–Rx separation can a polarimetric cluster tolerate before its polarimetric signatures
break down? This section documents the post-processing chain for that question. It covers what the
exported data contain, how the target return is isolated, which signatures are computed, and how
to interpret them. The verified routines are the module-level functions of
`polarimetry_processor.py` (pure NumPy, no AEDT), chained by its `MimoPolarimetryProcessor` class
and covered by `tests/test_polarimetry_processor.py`.

### 4.1 Setup

* **Source:** `DualPolHornCluster` (HFSS). Four open-ended WR-12 waveguides (TxH, TxV, RxH, RxV),
  each with its own CS. Wave ports sit 20 mm behind the apertures. λ₀ = 3.795 mm (79 GHz).
  *Phase Center Per Port* is set to each element's CS, so in SBR+ every port is a separate point
  source at its own position, with its own embedded far-field pattern.
* **Geometry:** the elements sit on a square of diagonal `copolarSpacing = copolarSpacingLambda·λ₀`
  around the origin, at (±d/2√2, ±d/2√2):

  | Element | Position (units of d/2√2) |
  | --- | --- |
  | TxV | (+1, +1) |
  | RxV | (−1, −1) |
  | TxH | (+1, −1) |
  | RxH | (−1, +1) |

  Co-polar pairs (TxH–RxH, TxV–RxV) are d apart along the two diagonals. Cross-polar pairs are
  d/√2 apart along x.
* **Scene:** the hand-built `LinkedDualPolHornCluster` (SBR+), which predates `sbr_simulator`. It
  contains the far-field linked antenna and a PEC sphere of radius `targetRadius` = 100 mm,
  centred at `targetRange` = 3 m on boresight. For designs built by `sbr_simulator --target sphere`,
  pass `targetDistance` as `--target-range` and `targetSize` as `--target-radius`. The setup
  uses 4 rays/λ and 5 bounces, and sweeps 74–84 GHz in 201 points (Δf = 50 MHz). The parametric
  sweep `copolarSpacingLambda` = 1.41…8.41 (step 1) gives a co-polar bistatic angle
  β = 2·atan(d/2R) of only 0.11°…0.63°.

### 4.2 Running it

```bash
# 1. Export the scattered-field solution (in the container via the `ansys-vnc link` wrapper,
#    attached to the GUI session; refuses to start an AEDT)
sbr_results mimo_polarimetry \
    LinkedDualPolHornCluster data/interim/bistatic_sweep/sphere_boresight_3m_scattered.npz \
    --setup "Setup : Sweep_Scattered" --sweep-variable copolarSpacingLambda
# 2. Analyse (host or container); figures + CSVs go to --out
bistaticity_analysis \
    data/interim/bistatic_sweep/sphere_boresight_3m_scattered.npz \
    --out data/processed/bistaticity_scattered
```

A solved SBR+ design holds three solutions under the same setup and sweep:

| Solution | Contents |
| --- | --- |
| `Setup : Sweep` | Total fields |
| `Setup : Sweep_Incident` | Incident only: the direct Tx→Rx coupling (§4.3) |
| `Setup : Sweep_Scattered` | Scattered only: the target |

Total = Incident + Scattered holds to 3·10⁻¹⁰ in the horn-cluster sweep. **Use
`Sweep_Scattered`.** The total solution also works in two ways:

* with `--background <Sweep_Incident or no-target run>.npz` for coherent subtraction, which
  reproduces the scattered-solution signatures to 5·10⁻⁹;
* via range gating (§4.4).

`--gate auto` (the default) decides between these. It gates only total-field exports without a
background. Target-only data (`_Scattered` exports, or a subtracted background) are processed
ungated over the full sweep, because gating would only discard band and smooth the response.
`--gate on|off` overrides the choice.

The same chain from Python:

```python
from pomaar.simulator.polarimetry_processor import MimoPolarimetryProcessor, REFERENCE_SPHERE

processor = MimoPolarimetryProcessor.from_npz("sphere_boresight_3m_scattered.npz")
processor.compensate_geometry(positions, focal_point)  # optional; per-variation element positions
signatures = processor.signatures(REFERENCE_SPHERE)    # dict of (n_variations, n_frequencies) arrays
```

**Quads.** The processor works on one polarimetric quad (RxH, RxV, TxH, TxV), found from the port
names as in `sbr_simulator`'s reports. Hand-built names (`HornCluster1_RxH`) and builder names
(`Linked1_Port_Rx_1_H`) both work. An export with several quads needs a choice:
`MimoPolarimetryProcessor.available_quads(path)` lists them (e.g. `Rx_1/Tx`, `Rx_2/Tx`), and
`from_npz(path, quad="Rx_2/Tx")` picks one. A dict `{"RxH": ..., "RxV": ..., "TxH": ..., "TxV": ...}`
of excitation names picks any four ports, e.g. single-pol Tx elements with different labels.

**Arrays from `hfss_array_builder`.** `sbr_results` attaches the builder's port table to the export
(below). `from_npz` then takes the port field signs out (`correct_signs=False` keeps the channels as
measured), so that every channel refers to +x (H) and +y (V); a sphere calibration could not reveal
a wrong sign. `compensate_geometry(None, focal_point)` uses the ports' phase centres from the
table. They are the values at build time, so a sweep of a layout variable still needs explicit
positions per variation.

`MimoPolarimetryProcessor` provides the following methods:

| Method | What it does |
| --- | --- |
| `gate` | Range gating (§4.4) |
| `calibrate_with_sphere` | Complex channel-imbalance calibration from a sphere measurement |
| `inject_noise` | Adds complex white Gaussian noise at a given SNR |

`calibrate_with_sphere` can only assume equal cross-pol channel gains, because a sphere does not
reveal them; see its docstring.

Tests (pure NumPy, no AEDT): `uv run --with pytest python -m pytest tests`.

`sbr_results.export_sbr_parametric_sweep` stores the following in the `.npz`:

* `s` (variation × channel × frequency, complex) and `frequencies_hz`;
* `sweep_values`;
* `expressions`, i.e. the AEDT names `S(<Rx port>,<Tx port>)`;
* the nominal design variables;
* for an array built by `hfss_array_builder`, its port table (`metadata["ports"]`): per excitation
  the builder port, role, polarization, field sign and phase centre (mm). `sbr_results` finds the
  `<ArrayDesign>_ports.json` that matches the excitations next to the project, or takes
  `--ports <file>`.

`sbr_results --port` attaches to a GUI or non-graphical AEDT session, but never starts one.

> [!NOTE]
> If `get_solution_data` returns *"No Data Available"*, the SBR+ results were invalidated (e.g.
> by an edit of the linked source design). Re-solve the SBR+ design; the exporter then works.

**Convention.** Scattering matrices are indexed `S[..., rx, tx]` with 0 = H and 1 = V, so
`S_HV` = receive H, transmit V = AEDT `S(...RxH, ...TxV)`.

### 4.3 What the exported S-parameters actually contain

The total-field SBR+ S-parameters (`Setup : Sweep`, the GUI default) are **not** the target
response. Transforming each channel to the two-way-path (delay) domain (`delay_profile`) shows two
components:

| Component | Two-way path (wrapped) | Level |
| --- | --- | --- |
| Direct Tx→Rx coupling between the proxy ports | d + 2 × 20 mm feeds ≈ 0.055–0.085 m | −35 dB (1.41λ) to −67 dB (8.41λ) |
| Sphere specular return | 2(R − a) + feeds = 5.85 m, wraps to −0.145 m | ≈ −97.5 dB |

* **The wrap.** The unambiguous path of the sweep is c/Δf = 5.996 m, so the sphere sits only
  ≈ 0.2 m from the coupling peak. That is about 7 resolution cells of c/B = 3 cm.
* **The level check.** −97.5 dB matches the radar equation for σ = πa² at 2.9 m with about
  7.5 dBi element gain. Plotting the raw S-parameters therefore plots the coupling, 30–60 dB above
  the target.
* **Dispersive coupling.** The coupling is not a clean point response. At 3.41λ and 5.41λ its
  phase swings by about 80° across the band.
* **Periodic kinks.** All variations show small kinks every ≈ 1.25 GHz (75.25, 76.5, 77.8, 79.0,
  80.3, … GHz). These produce paired echoes 0.24 m (c / 1.25 GHz) from every peak, about 31 dB
  down. They are common to all channels, so ratios are unaffected.

Because the coupling is 30–60 dB stronger, a 0.1 % irregularity in it reaches the target's level
at every delay. When the total solution is gated, this contaminates the co-pol channels at 3.41λ
and 5.41λ. It also leaks into the cross-pol channels *inside* the gate, where the guard-band check
cannot see it: gated total-field cross-pol came out at −30…−41 dB at 79 GHz, whereas the scattered
solution gives −36…−57 dB (§4.7).

> [!IMPORTANT]
> **Remove the coupling at its source.** Export `Setup : Sweep_Scattered`. Equivalently, subtract
> `Sweep_Incident` or a no-target run (`Sphere1` → *Model* off, or `sbr_simulator --target none`)
> with `--background`. Gating the total field is only a fallback, and its cross-pol values are not
> trustworthy.

### 4.4 Total-field fallback: gate-normalized time gating

Everything in this subsection applies only when gating is on (total-field data without a
background).

`time_gate(s, f, centre_path, half_width, window)` works as follows:

1. Taper the sweep with `window`.
2. Transform to the path domain.
3. Keep a rectangular gate of ±`half_width` around the (wrapped) target path.
4. Transform back.
5. Divide by the same processing applied to an ideal point reflector at `centre_path`, then
   multiply that reflector back in (gate normalization).

Step 5 removes the taper and the truncation of the target's own main lobe. What remains is only
leakage from other paths.

**Chosen parameters** (validated on synthetic data with the measured levels and paths):

| Parameter | Value | Why |
| --- | --- | --- |
| `window` | Kaiser β = 14 | Hann leaks the −35 dB coupling into the gate (6 dB / 43° error) |
| Half-width | ±0.06 m | Narrower truncates the target; wider reaches the coupling |
| Valid band | taper ≥ 0.5 → 77.45–80.55 GHz (`gate_valid_band`) | Errors grow outside it |

**Synthetic errors:**

| Case | Error |
| --- | --- |
| Sphere, worst case (−35.6 dB coupling) | ≤ 0.015 dB / 0.12° |
| −140 dB cross-pol term at the same delay | ≤ 1 dB / 7° |

**Gate centre.** The gate is centred on the peak of |S_HH|² + |S_VV|² in the zero-padded profile,
searched within ±3 cm of the geometric path. It is common to all four channels, so ratios see the
same gate.

**Contamination checks** (`delay_window_peak`). Each check marks points as hollow markers in the
figures:

* **Gate SNR:** the target peak in the gate versus the peak in a guard band at −0.40…−0.18 m from
  it, on the side away from the coupling. Clean spacings reach about 31 dB, limited by the paired
  echo. Co-pol below 20 dB is flagged.
* **Peak shift:** a co-pol peak more than 5 mm from the sweep's median target path is flagged
  (3.41λ: +17 mm, 5.41λ: +15 mm).
* **Cross-pol floor:** the guard-band level of the cross-pol channel relative to the co-pol
  target, drawn dotted in panel (3). Cross-pol less than 10 dB above it is "unresolved", and every
  metric using cross-pol channels is flagged.

### 4.5 Signatures (panels of `signatures_vs_spacing.png`)

Before the signatures are computed, `compensate_geometry` divides each channel by its
spherical-wave phase to the sphere's front point (`spherical_wave_phase`), gated or not. For this
square every element is d/2 from the centroid, so all four channels have identical paths (spread
0 µm), and compensation changes no ratio (≤ 3·10⁻¹⁴ °). It matters for off-boresight targets
(§4.6).

Each panel plots one scalar against spacing, with β on the top axis. The line is the 79 GHz value
and the band is the min–max over the valid band: the whole sweep when ungated, 77.45–80.55 GHz
when gated.

| # | Signature | Function | Sphere ideal | Dihedral (edge ∥ H) ideal |
| --- | --- | --- | --- | --- |
| 1 | Co-pol imbalance \|S_VV/S_HH\| | `copolar_ratio` | 0 dB | 0 dB |
| 2 | Co-pol phase difference ∠(S_VV S_HH*) | `copolar_ratio` | 0° (reference) | 180° relative to sphere |
| 3 | Cross-pol leakage per Tx pol: \|S_VH/S_HH\| (Tx H), \|S_HV/S_VV\| (Tx V) | `crosspolar_isolation` | → −∞ | → −∞ |
| 4a/b | Reciprocity S_HV/S_VH | `reciprocity_ratio` | (undefined, both → 0) | (use a 45° dihedral) |
| 4c | Normalized asymmetry \|S_HV − S_VH\| / ‖S‖_F | `normalized_asymmetry` | → −∞ | → −∞ |
| 5 | Bistatic Pauli power fractions of k = [HH+VV, HH−VV, HV+VH, j(HV−VH)]/√2 | `pauli_vector`, `pauli_power_fractions` | k₁ = 1 | k₂ = 1 |
| 6 | Infidelity 1 − \|⟨k, k_ref⟩\| / (‖k‖‖k_ref‖) | `polarimetric_fidelity` | → −∞ | → −∞ (vs `REFERENCE_DIHEDRAL`) |
| 7 | Single-target Cloude α, cos α = \|k₁\|/‖k‖; entropy H of the band-averaged T | `alpha_angle`, `coherency_matrix`, `entropy_alpha` | 0°, H = 0 | 90°, H = 0 |

Notes on the table:

* **Antisymmetric component.** k₄ is zero for any monostatic reciprocal measurement. Its power
  fraction is the direct bistaticity indicator. Note that the monostatic k_p of the literature
  (`bistatic=False`) symmetrizes it away.
* **Entropy.** H is computed for the 4×4 (bistatic) and 3×3 (symmetrized) coherency matrices
  averaged over the valid band. For one deterministic target, H > 0 means its signature changes
  across the band.
* **Phase reference.** All phases are relative between channels. The absolute arg S_HH is path
  length (geometry), not polarimetry.

### 4.6 The linear term versus what survives compensating it

For a Tx at p_t and Rx at p_r, a far-field target in direction r̂ adds the phase
**k (p_t + p_r) · r̂**. This is the steering phase of the virtual element at (p_t + p_r)/2. In
this cluster the four channels have different virtual element positions:

| Channel | Pair | Virtual element |
| --- | --- | --- |
| HH | TxH + RxH | (0, 0) |
| VV | TxV + RxV | (0, 0) |
| HV | TxV + RxH | (0, +d/2√2) |
| VH | TxH + RxV | (0, −d/2√2) |

So HH and VV share a phase centre, while HV and VH are displaced along y in opposite directions.
For a target at elevation θ_y off boresight (y–z plane):

**Δφ(HV − VH) = √2 · k · d · sin θ_y**, and **Δφ(HV − HH) = k · d/√2 · sin θ_y**.

This breaks the apparent reciprocity S_HV = S_VH. Take a 45°-rotated dihedral (S_HV = S_VH,
HH = VV = 0) seen through the uncompensated linear term. A fraction sin²(Δφ/2) of its power moves
from k₃ to the antisymmetric k₄, which helix-sensitive decompositions (Yamaguchi, Singh) read as
helix scattering. Values from `linear_term.csv` and the first two panels of
`linear_term_study.png`:

| d | θ_y | Δφ(HV − VH) | 45° dihedral power into k₄ |
| --- | --- | --- | --- |
| 1.41λ | 1° | 12.5° | 1.2 % |
| 1.41λ | 2° | 25.1° | 4.7 % |
| 1.41λ | 5° | 62.6° | 27 % |
| 1.41λ | 10° | 124.7° | 78 % |
| 4.41λ | 1° | 39.2° | 11 % |
| 4.41λ | 2° | 78.4° | 40 % |
| 8.41λ | 1° | 74.7° | 37 % |
| 8.41λ | 2° | 149.4° | 93 % |

The k₄ fraction is periodic in θ_y, so it can also vanish by coincidence.

**The linear term is deterministic and removable.** Divide each channel by `plane_wave_phase`
for its own Tx/Rx pair and the look direction. Delay-and-sum beamforming does this implicitly, as
long as the steering vectors use each polarization channel's own virtual positions, not one shared
array. Better still, divide by `spherical_wave_phase` for the focal point, which also removes the
Fresnel (quadratic) term. **Polarimetric errors caused by the linear term are therefore a
processing error, not a limit of the configuration.**

**What survives compensation** sets the real bistaticity limit:

1. **Fresnel residual, when only plane-wave compensation is used** (third panel of
   `linear_term_study.png`). Elements on a circle around the centroid make it common-mode at
   boresight. For d = 8.41λ the differential residual HV − VH stays ≤ 0.95° at R = 0.5 m and
   θ_y ≤ 30°, ≤ 0.24° at 1 m, and ≤ 0.03° at 3 m. Negligible here, and zero with spherical
   compensation.
2. **Aspect-angle diversity.** Each Tx/Rx pair sees the target under its own bistatic geometry, so
   an extended target presents a genuinely different S per channel. The relevant parameter is β
   compared with the target's angular lobe width λ/L, i.e. d·L/(λR). A sphere is insensitive
   (specular at any β); flat and retro-reflecting targets are sensitive (a dihedral with
   L ≈ 100 mm has λ/L ≈ 2.2°, a trihedral more so).
3. **Antenna diversity.** HV uses TxV + RxH while VH uses TxH + RxV, so the two cross-pol channels
   are measured by different antenna pairs. Their embedded patterns differ in gain, phase and
   cross-pol (mutual coupling to the orthogonal neighbour). The off-boresight angle also differs
   per element. None of this cancels under geometric compensation. It needs polarimetric
   calibration (per-channel distortion matrices from a sphere and a depolarizing reference).
4. **Coupling and multipath.** Direct coupling is removable (scattered solution, background
   subtraction or gating).
   Interactions between the target and the antenna structure are absent from this SBR+ model
   unless the structure is included in the scene.

### 4.7 First results: boresight sphere at 3 m

From `Setup : Sweep_Scattered`, processed ungated (default), i.e.
`data/processed/bistaticity_scattered/signatures_vs_spacing.csv` (79 GHz values):

| d (λ₀) | β (°) | \|S_VV/S_HH\| (dB) | ∠ (°) | Tx H leak (dB) | Tx V leak (dB) | \|S_HV/S_VH\| (dB) | Infidelity (dB) | α (°) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.41 | 0.11 | −0.08 | −0.4 | −41.3 | −36.5 | 4.7 | −40.4 | 0.78 |
| 2.41 | 0.18 | −0.07 | 0.8 | −40.6 | −57.5 | −16.9 | −42.4 | 0.62 |
| 3.41 | 0.26 | −0.14 | 0.6 | −40.3 | −40.5 | −0.4 | −40.4 | 0.77 |
| 4.41 | 0.33 | −0.19 | 0.8 | −41.1 | −48.5 | −7.6 | −39.7 | 0.84 |
| 5.41 | 0.41 | 0.00 | 0.6 | −42.9 | −44.2 | −1.2 | −44.5 | 0.48 |
| 6.41 | 0.48 | −0.12 | −0.8 | −39.1 | −47.9 | −8.9 | −40.7 | 0.74 |
| 7.41 | 0.56 | 0.15 | 0.1 | −43.5 | −41.4 | 2.2 | −41.5 | 0.68 |
| 8.41 | 0.63 | −0.02 | −0.1 | −44.2 | −42.8 | 1.4 | −46.3 | 0.39 |

* **Co-pol.** Ideal: within −0.28…+0.20 dB and −1.3…+2.7° over the full 74–84 GHz band, with no
  bistatic trend up to β = 0.63°, as expected for a sphere. Absolute |S_HH| is −97.6…−98.1 dB.
* **Cross-pol.** Leakage is −33…−73 dB over the band (−36…−57 dB at 79 GHz), with no trend in β.
  It varies non-monotonically from one spacing to the next, which points to the elements'
  embedded-pattern cross-pol (mutual coupling to the orthogonal neighbour changes with spacing;
  item 3 of §4.6) rather than to the bistatic geometry.
* **HV versus VH.** They differ by up to 17 dB in magnitude and by arbitrary phase. With cross-pol
  produced by different antenna pairs, nothing forces them to agree. The antisymmetric k₄
  fraction (−40…−52 dB at 79 GHz) therefore reflects antenna diversity, not bistaticity.
* **Overall.** Infidelity ≤ −34 dB over the band, α ≤ 1.5° and H < 10⁻³ are all set by that
  cross-pol.
* **Interim conclusion.** Up to 8.41λ at 3 m (β ≤ 0.63°), the boresight sphere shows no
  polarimetric breakdown. The configuration's polarimetric floor (about −40 dB infidelity) is set
  by the antennas. To find where bistaticity starts to matter, use targets with narrow bistatic
  lobes, shorter range, or off-boresight positions (§4.8).
* **Gating versus ungated.** Gating the scattered data changes co-pol by ≤ 0.03 dB / 0.2° at
  79 GHz. Cross-pol changes by up to 4.5 dB, because the gate smooths over a few GHz, including
  the ≈ 1.25 GHz kinks of §4.3.
* **Total-field results.** The gated total-field analysis (`data/processed/bistaticity_total/`)
  flags 3.41λ and 5.41λ as contaminated. It agrees on co-pol at the clean spacings, but overstates
  the cross-pol (−30…−41 dB) because coupling leaks inside the gate.

### 4.8 Next steps

* **Finer frequency steps.** Halving Δf, which the SBR+ sweep inherits from the HFSS
  interpolating sweep, doubles the unambiguous path. It is only needed when gating total fields.
* **Cross-pol cross-check.** Check the cross-pol level against the HFSS embedded patterns:
  M_VH ≈ e_rx,V · S_ideal · e_tx,H at θ ≈ β/2 for every element.
* **Off-boresight targets.** Put the target off boresight (θ_y scan) to observe the linear term in
  simulation, and verify that per-channel compensation removes it.
* **Dihedral and trihedral.** `sbr_simulator --target dihedral|trihedral` builds both; use
  `targetRoll` for the edge parallel and perpendicular to each baseline and for the 45° dihedral
  (reciprocity), and a range R ≈ 1 m so that d·L/(λR) reaches order 1.
* **Generalize `bistaticity_analysis`.** It is still specific to this study: the sphere reference
  (`REFERENCE_SPHERE`), a boresight focal point at the sphere's front, and the `DualPolHornCluster`
  geometry (`cluster_positions`, 20 mm feeds). Dihedral, trihedral and off-boresight runs need the
  reference matrix, the target position and the element positions as inputs (e.g. read from the
  SBR+ design variables and the layout).
