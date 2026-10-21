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
ansys-vnc exec hfss_array_builder <path_to_project.aedt> <layout.yaml>
```

To run the scripts from any directory, link them once: `ansys-vnc link --host bistaticity_analysis ~/Repositories/pomaar` writes wrappers into `~/.local/bin`. `hfss_array_builder`, `sbr_simulator` and `sbr_results` then run in the container, and `bistaticity_analysis` runs on the host. Rerun the command after adding console scripts.
`ansys-vnc exec` keeps the current directory, pins the process to the compute cores and puts the nearest `.venv` (here the project's) first on `PATH`, so the console scripts and `python` run by name. `hfss_array_builder`, `sbr_simulator` and `sbr_results` all attach to the running AEDT GUI session (its gRPC port is found automatically) and never start an AEDT of their own: start the GUI first (`ansys-vnc aedt`) and open the project. The builder and `sbr_simulator` work in that session in plain view and save the project.

### Option B: Running natively (Windows / supported Linux)

Ensure that:

1. `ANSYSEM_ROOT252` (or your corresponding AEDT installation path) is set in your environment variables.
2. The `pomaar` package dependencies are installed (`uv sync`).
3. On Linux outside a container, set `POMAAR_ALLOW_NATIVE_AEDT=1` to let the scripts connect to the running AEDT GUI session.

---

## 2. Running the pipeline

### Module 1: HFSS full-wave array synthesis (`hfss_array_builder.py`)

This script automates the creation of a planar MIMO array on a single contiguous PCB board by copying, replicating, and boolean-cutting template geometries from an isolated unit-cell element design.

To run (prefix with `ansys-vnc exec` in the container setup):

```bash
hfss_array_builder <path_to_project.aedt> <layout.yaml> [--output-design NAME] [--use-existing-cs] [--overwrite] [--simulate | --no-simulate]
```

* **Example:**

  ```bash
  hfss_array_builder ~/Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt layouts/ApertureCoupledPatch_polarimetric_cluster.yaml
  ```

* **Inputs:** the project and a layout file; nothing else describes the array. A layout is bound to the unit-cell designs it is built from (the generators write them, see Layouts below). It holds:

  | Key | Meaning |
  | --- | --- |
  | `metadata.primary_unit_cell` (required) | The design every element uses unless told otherwise. It provides the analysis setup and names the array design `<primary_unit_cell><LayoutFileStem>`, where a file stem starting with `<primary_unit_cell>_` drops that prefix (`ApertureCoupledPatch_polarimetric_cluster.yaml` → `ApertureCoupledPatchPolarimetricCluster`). |
  | `metadata.unit_cells` | Further unit cells per polarization set: `h`/`v` for single-pol elements, `hv` for dual-polarized ones. |
  | element `unit_cell` | One element's own unit cell. |
  | `metadata.drop_unused_ports` | `true` removes the unused ports of dual-port cells (see below). |
  | `metadata.phase_reference_mode` | `per_port` (default) or `single`; see "What it does". |
  | `metadata.variables`, `metadata.board` | Design variables of the layout (spacings, `pcbMargin`, …) and the board size formulas in them; the builder adds `unitCellExtentX/Y`, `unitCellExtent` (the unit cells' reach from their phase centre) and `airboxClearance`. |
  | element `position`, `position_expression` | The numeric position (mm, x y z) drives the builder's checks; the parametric `[x, y]` expression places the element in AEDT. |

  Older layouts' `centre_frequency_ghz`, `bandwidth_ghz` and `topology` are ignored (with a note).

* **Options**, as in `sbr_simulator`: `--output-design` names the array design (default `<primary_unit_cell><LayoutFileStem>`); the others only answer the interactive prompts. `--use-existing-cs` reuses the unit cells' `PhaseCentreCS` instead of offering a fit, `--overwrite` replaces an existing array design (asked otherwise; never without a terminal), and `--simulate` / `--no-simulate` decide silently whether to solve the setup after building (asked otherwise; not solved without a terminal). The far-field spheres and reports are always set up.
* **Frequencies** are never part of a layout. The builder takes the centre frequency from the primary unit cell's setup (for the far-field reports, the λ/4 Airbox clearance and a phase-centre fit) and copies that setup as `Setup` for its meshing: adaptive frequencies, convergence and solver settings, as validated on the unit cell (a native copy-paste, else a recreation from its properties, else the build fails). It then replaces the copied sweeps by its own **discrete** sweep `Setup : Sweep` over the frequencies of the unit cell's (first) sweep, saving the radiated fields at every frequency, and `sbr_simulator` links that sweep and copies its frequencies into the SBR+ setup. The sweep is discrete because SBR+ interpolates a linked antenna's far fields between the frequencies where the HFSS design saved fields (Module 2); the price is one full-wave solve of the array per frequency, on the final adaptive mesh. The unit cell itself can keep any sweep type; it needs a sweep, and the builder warns when its adaptive frequencies stay below the top of that sweep. `ApertureCoupledPatch` sweeps 77–81 GHz in 81 points, i.e. Δf = 50 MHz: the unambiguous two-way path is c/Δf = 6.0 m, a maximum unambiguous range of c/(2Δf) = 3.0 m, and the range resolution is c/(2B) = 3.7 cm. A target beyond 3 m (e.g. `sbr_simulator`'s default 4 m) folds back in the range profiles and maps; the signatures, taken in the frequency domain, are unaffected. Change the sampling in the unit cell (e.g. 161 points for 6 m).
* **Polarization:** an element's polarization comes from its unit-cell template, not from its placement. Port sheets carry their native polarization in the name (`PortSheet_H`, `PortSheet_V`); a unit cell has one such port (single-pol) or two (dual-port). Elements without a unit cell of their own polarization set use the primary design, turned by their `yaw`. Either way the builder checks that R(yaw) applied to each port's native E-field lies along x for `polarization: h` and y for `v` (`hv` needs a dual-port cell) and aborts otherwise, before it touches any design.
* **Unused ports:** an element that uses only one port of a dual-port cell keeps the other feed, terminated in a matched wave port `Port_<label>_<P>_Load`, so that it does not become an open stub; `metadata.drop_unused_ports: true` removes it instead. `sbr_simulator` skips `_Load` ports.
* **Field signs:** turning an element by 180° flips its E-field. After building, the builder prints every port's polarization and sign (+1 along +x for H, +y for V) and writes them, with the phase centres, to `<ArrayDesign>_ports.json` next to the project. The measured matrix is diag(s_RxH, s_RxV) S diag(s_TxH, s_TxV), so a role whose H and V ports have opposite signs (pair sign −1) negates the cross-pol channels against the co-pol ones; a sphere calibration cannot reveal this, and the builder warns about it.
* **What it does:** creates or updates the array design, syncs design variables, sizes the PCB board to fit the feedlines, replicates the elements, assigns wave ports (lumped-port fallback), excites the first Tx port for far-field reports, draws a $\lambda_0/4$ vacuum radiation Airbox, checks every element's placement and that every wave port lies on the Airbox, and sets up Theta-Phi far-field spheres in `PhaseCentreCS` with realized-gain pattern and directivity-based XPD reports. The far-field phase reference is set per port to its own CS (`phase_reference_mode: per_port`), so that a linked SBR+ antenna launches every port's rays from its own phase centre; `single` references all ports to the array centroid `PhaseCentreCS` instead, which is only valid for targets in the far field of the whole array.

#### Layouts (`layouts/`)

Each YAML is written by the `generate_<topology>.py` next to it (`--help` lists the options) and
named after its primary unit cell, `<UnitCell>_<topology>.yaml`. Every generator takes:

| Option | Layout metadata |
| --- | --- |
| `--unit-cell <Design>` (required) | `primary_unit_cell` |
| `--extra-unit-cell h\|v\|hv=<Design>` (repeatable) | `unit_cells` |
| `--drop-unused-ports` | `drop_unused_ports: true` |
| `--phase-reference-mode per_port\|single` | `phase_reference_mode` |
| `--centre-freq <GHz>` | none: only the design frequency of the λ-based default spacings |
| `-o <file>` | another output file; its stem then names the array design |

So the same topology for several unit cells is one generator run per cell, each followed by the
builder on the file it names:

```bash
for cell in ApertureCoupledPatch SeriesFedPatch; do
    python layouts/generate_polarimetric_cluster.py --unit-cell $cell
    ansys-vnc exec hfss_array_builder <project.aedt> layouts/${cell}_polarimetric_cluster.yaml --no-simulate
done
```

The shipped YAMLs are bound to `ApertureCoupledPatch`, a V-polarized cell whose feed points to −y
at yaw 0, and regenerate exactly with:

```bash
python layouts/generate_contiguous_ula.py --unit-cell ApertureCoupledPatch --tx-offset-y 10               # _contiguous_ula_4x4
python layouts/generate_contiguous_ula.py --unit-cell ApertureCoupledPatch --tx-count 3 --tx-offset-y 10  # _contiguous_ula_3x4
python layouts/generate_spatial_diversity_2d.py --unit-cell ApertureCoupledPatch
python layouts/generate_polarimetric_cluster.py --unit-cell ApertureCoupledPatch
python layouts/generate_polarimetric_ula.py --unit-cell ApertureCoupledPatch [--extra-unit-cell hv=<DualCell>]
```

| Layout | Elements | Status (with `ApertureCoupledPatch`) |
| --- | --- | --- |
| `…_contiguous_ula_3x4`, `…_contiguous_ula_4x4` | V-only Rx and Tx ULAs | Builds; single-pol, so no polarimetric quad |
| `…_polarimetric_cluster` | Tx/Rx × H/V, H by turning the V cell (feeds pointing outwards) | Builds; pair signs +1, one quad `Rx/Tx` |
| `…_spatial_diversity_2d` | V-only L-shaped Rx and triangular Tx | Does not build: the elevation Rx and the outer Tx sit behind other elements, so their feeds end inside the board, where a wave port cannot be (the placement check aborts) |
| `…_polarimetric_ula` | Dual-polarized `hv` Rx elements, alternating V/H Tx | Needs a dual-port cell (`--extra-unit-cell hv=<Design>`); two Tx per polarization, so the processor needs the four ports of a quad explicitly |

The earlier `polarimetric_ula` placed a V and an H copy of the single-pol cell at every Rx position. That geometry overlaps and gave both roles a pair sign of −1, so it is no longer generated.

### Module 2: SBR+ target solver (`sbr_simulator.py`)

This script builds an SBR+ design next to a solved HFSS antenna design of the same project: it places a canonical PEC target, links the HFSS design as the antenna source, pairs every Tx port with every Rx port, copies the frequency sweep, creates quick-look S-parameter reports and optionally solves.

To run (prefix with `ansys-vnc exec` in the container setup; `uv sync` installs the script):

```bash
sbr_simulator <path_to_project.aedt> <hfss_design_name> [--boresight x|z] [--target dihedral] [--size 40mm] [--range 4m] [--azimuth 0deg] [--elevation 0deg] [--roll 0deg] [--map-variables VAR ...] [--field-type farfield] [--simulate | --no-simulate]
sbr_simulator <path_to_project.aedt> <sbr_design_name> --scene-only [--target ...] [--range ...] [--overwrite]
```

* **Example:**

  ```bash
  sbr_simulator ~/Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt DualPolHornCluster --target sphere --map-variables copolarSpacingLambda
  ```

* **Units:** `--range` (alias `--distance`) is in metres and `--size` in mm unless a unit is given (`4m`, `4 metres`, `4 meters` and `400mm` all work); metres are written to AEDT as its unit token `meter`.
* **What it does**, in the design `Linked<hfss_design_name>` (`--output-design` overrides; an existing one is only replaced after confirmation or `--overwrite`):
  1. Creates `AntennaCS` at the global origin and places the linked antenna in it. With the default `--boresight x` (ISO 8855: x forward, y left, z up), `AntennaCS` is turned by Euler angles ZYZ (0°, 90°, 90°), so the antenna's own frame (boresight +z, H along x, V along y) radiates along +x with H along +y and V along +z; with `--boresight z` it is the global CS.
  2. Creates `TargetCS` at range `targetRange`, azimuth `targetAzimuth` and elevation `targetElevation`. With x boresight it is defined in the global CS at `X = targetRange cos(targetElevation) cos(targetAzimuth)`, `Y = targetRange cos(targetElevation) sin(targetAzimuth)`, `Z = targetRange sin(targetElevation)` (azimuth from +x towards +y, elevation above the xy plane), with Euler angles ZYZ `(Phi, Theta, Psi) = (targetAzimuth, 90deg - targetElevation, 90deg + targetRoll)`: its local z axis points radially away from the antenna and, at zero roll, its local x axis lies along H and its local y axis along V at any azimuth and elevation, so `targetRoll` is a pure rotation about the line of sight (e.g. a 45° dihedral). With z boresight (the earlier convention) it sits in `AntennaCS` at `X = targetRange cos(targetAzimuth) sin(targetElevation)`, `Y = targetRange sin(targetAzimuth) sin(targetElevation)`, `Z = targetRange cos(targetElevation)`, `targetElevation` being the polar angle off boresight, with Euler angles `(targetAzimuth, targetElevation, targetRoll - targetAzimuth)`. Edit the TargetCS angles directly for other orientations.
  3. Creates the target in `TargetCS`, facing the antenna along local -z, sized by `targetSize` (`--size`) and made PEC by a Perfect E boundary (sheets carry no material); its SBR+ surface roughness and height standard deviation are 0 unless `--roughness` / `--height-deviation` are given:
     * sphere: PEC solid of radius `targetSize`, centred at `TargetCS`;
     * dihedral: two square `targetSize` × `targetSize` sheets, i.e. `targetSize` is the seam length and the sheet depth. The seam lies on the local x axis through `TargetCS` (parallel to H at zero roll), which is also the phase reference of the double bounce. One sheet is drawn, tilted by 45° and copied by a 180° rotation about the symmetry axis (local z);
     * trihedral: three right-isosceles triangular sheets with edge length `targetSize`, the corner at `TargetCS` and the symmetry axis pointing at the antenna. One sheet is drawn and copied by 120° rotations about local z.

     The defaults (40 mm at 4 m) keep all three in their own far field at 79 GHz (2D²/λ ≤ 3.4 m); the script warns when the chosen range is shorter. `--target none` builds the same scene without a target; its solve is the background (the direct Tx→Rx coupling between the proxy ports), a cross-check of the incident solution below (use `--output-design` to keep both designs).
  4. Links the HFSS design (solution `--source-solution`, default its first discrete sweep that saves fields, e.g. the builder's `Setup : Sweep`; anything else is refused before an existing SBR+ design is touched, see below) with *Simulate source design as needed*, *Model Visualization* of all its objects except the radiation boundary objects (`RadiatingSurface`, the builder's `Airbox`), and far-field (standard) or near-field (bumper/radome integration) results, asked interactively unless `--field-type` is given. `--map-variables` maps independent HFSS design variables to SBR+ design variables of the same name, which are created from the HFSS value when missing, so they can be swept in the SBR+ design. Each mapped variable is also declared as a parameter of the linked component's definition, as the GUI's linked-antenna dialog does; without that declaration (links made by `sbr_simulator` before this fix), AEDT keys the link's per-instance caches by no parameter at all: changing the variable resizes the footprint box but not the visualized antenna, and a horn link made that way returned the same S-parameters at every swept value. Re-create such links (GUI or `sbr_simulator`).
  5. Selects Tx/Rx: ports whose names contain `Tx`/`Rx` and end in `H`/`V` (e.g. `TxH`, `Port_Tx_1_V`); every Tx transmits to every Rx. The builder's matched `_Load` ports are skipped.
  6. Creates `Setup` (4 rays per wavelength, 5 bounces, PTD + UTD, creeping waves) whose sweep `Sweep` copies every range of the linked HFSS sweep, so SBR+ runs exactly at the frequencies where the HFSS design saved fields.
  7. Creates quick-look reports from the scattered solution (see below): `S-parameters dB` for all Tx/Rx pairs and, per polarimetric quad, `Co-pol VV/HH` (`S_VV/S_HH` in dB and degrees; ideally 0 dB for both, with the dihedral's phase 180° away from the sphere's) and `Cross-pol isolation` (`S_VH/S_HH` for Tx H, `S_HV/S_VV` for Tx V, in dB). `S-parameters dB (total)` shows the total solution for reference. A quad pairs a dual-polarized Rx with a dual-polarized Tx: an element with an H and a V port of the same name prefix (`Port_Rx_1_H`/`_V`), or a role's single H and V ports when it has exactly one of each (the polarimetric cluster). When the HFSS design comes with the builder's port table (`<design>_ports.json` next to the project), the VV/HH phase is corrected for the port field signs.
  8. Asks whether to solve (`--simulate` / `--no-simulate` decide silently). Solving also solves the HFSS design where its linked solution is missing.

  `SbrSimulationManager.export_results(path, sweep_variable)` writes the solved scattered solution to an `.npz` for Module 3 (via `sbr_results.py`, which is also a standalone CLI; see §4.2).

* **Rebuilding the target scene of an existing SBR+ design** (`--scene-only`, `SbrSimulationManager.rebuild_target_scene`): for an SBR+ design that has its antennas already, e.g. parametric-beam sources, only `TargetCS` and the target are (re)built as above; the antennas, Tx/Rx pairs, setup and reports are kept. An existing canonical target (`TargetSphere`, `TargetDihedral`, `TargetTrihedral` with its `…PerfectE` boundary), `TargetCS` and the former range variable `targetDistance` are replaced after confirmation (or `--overwrite`).
* **Setup and sweep names are read back.** The setup and sweep are requested as `Setup` and `Sweep`, but AEDT may de-duplicate names (`Setup1`, …). The reports, the solve and `export_results` use the names AEDT actually gave them (`analysis_solution`): the setup created in this run, else the design's only (or first) setup and its first sweep. `sbr_results` likewise defaults to the design's own `<setup> : <sweep>_Scattered`. The examples below write `Setup : Sweep`.
* **Only discrete sweeps are linked.** SBR+ receives the linked antenna's far fields only at the frequencies where the HFSS design saved fields, and interpolates them in between. Linking `Setup : LastAdaptive` gives fields at the adaptive frequencies only; an interpolating sweep gives them at its basis frequencies only, which it picks by bisecting the band. Both put ripple into every SBR+ result, with slope breaks exactly at the solved frequencies: the dihedral export of §4 peaked at 77, 79 and 81 GHz (its link was `LastAdaptive`), and the sphere export at 75.25, 76.5, …, 82.75 GHz (an interpolating sweep), with curvature 100–250 times the median at those frequencies. `sbr_simulator` therefore links only a discrete sweep that saves fields (all fields or the radiated ones) and refuses anything else, naming each candidate and why. It reads the sweeps from the saved project file after saving the project, since PyAEDT 1.1 reports every loaded sweep as interpolating. A hand-built HFSS design needs a discrete sweep with saved fields over the SBR+ band.
* **The scattered field is the default, everywhere.** SBR+ splits its total solution `Setup : Sweep` into `Setup : Sweep_Incident`, the direct Tx→Rx coupling between the linked antenna's ports (the line-of-sight path that never touches the target), and `Setup : Sweep_Scattered`, everything that reached the receiver via the target (total = incident + scattered). The reports above, `export_results` and `sbr_results` all use the scattered part. This is the simulation's counterpart of what every practical radar measurement does to the antenna coupling: the direct Tx→Rx leakage arrives at (almost) zero range, well ahead of any target, and is removed by range gating, i.e. by keeping only the target's cells of the range (delay) profile, or equivalently by subtracting an empty-scene background measurement. SBR+ makes the same separation exactly, by bookkeeping of the rays instead of by delay, so it needs no window, loses no bandwidth or resolution, and still works when the frequency step folds the target's delay next to the coupling (§4.4). Because the coupling is 30–60 dB above the target at short Tx/Rx spacings, reading the total field means reading the coupling. The total field stays visible in `S-parameters dB (total)` and exportable with `sbr_results --setup "Setup : Sweep"` for debugging; the processor and the analysis assume the scattered solution.

### Module 3: Offline polarimetric post-processing (`polarimetry_processor.py`, `bistaticity_analysis.py`)

Pure NumPy, no AEDT. `MimoPolarimetryProcessor` runs the usual MIMO radar chain on an exported sweep (frequency window, IFFT to range profiles, optional spatial window, beamforming with steering vectors, maps of total power, polarimetric channel powers and Pauli features); `CanonicalTargetProcessor` adds the qualities of canonical scatterers; `bistaticity_analysis` plots those against the swept variable. Section 4 documents the chain, the signatures and the first results.

### Spelling

The code uses British spelling with Oxford -ize (centre, metre, polarization). Command-line flags also accept the American spelling (`--use-existing-phase-center`, `--center-freq`), and names read from AEDT designs are looked up in both spellings (e.g. a unit cell's `PhaseCentreCS` or `PhaseCenterCS`, `--map-variables centreWavelength` for an HFSS `centerWavelength`). Lengths read from AEDT accept `meter(s)` and `metre(s)`.

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
> 5. **Several unit cells in one array:** all of them must have the same board stack (layer names, materials and z extents) and must not define a design variable differently (rename such variables per unit cell). The analysis setup comes from the primary unit cell, and the board stack is placed at its phase-centre height.

---

## 4. Bistaticity study: polarimetric signatures from SBR+ sweeps

How much Tx–Rx separation can a polarimetric cluster tolerate before its polarimetric signatures
break down? This section documents the post-processing chain for that question: what the exported
data contain, how they are processed, which signatures are computed, and how to interpret them.
The routines are the functions of `polarimetry_processor.py` (pure NumPy, no AEDT), chained by its
`MimoPolarimetryProcessor` and `CanonicalTargetProcessor` classes and covered by
`tests/test_polarimetry_processor.py`.

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
  centred at `targetRange` = 3 m on boresight. The setup uses 4 rays/λ and 5 bounces, and sweeps 74–84 GHz in 201 points (Δf = 50 MHz). The parametric
  sweep `copolarSpacingLambda` = 1.41…8.41 (step 1) gives a co-polar bistatic angle
  β = 2·atan(d/2R) of only 0.11°…0.63°.

### 4.2 Running it

```bash
# 1. Export the scattered-field solution, the default (in the container via the `ansys-vnc link`
#    wrapper, attached to the GUI session; refuses to start an AEDT)
sbr_results mimo_polarimetry \
    LinkedDualPolHornCluster data/interim/bistaticity/horn_cluster_sphere_boresight_3m.npz \
    --sweep-variable copolarSpacingLambda
# 2. Analyse (host or container); figures + CSV go to
#    data/processed/bistaticity/horn_cluster_sphere_boresight_3m/
bistaticity_analysis data/interim/bistaticity/horn_cluster_sphere_boresight_3m.npz
```

`bistaticity_analysis` options: `--azimuth`/`--elevation` (deg; azimuth from boresight towards H,
elevation towards V) steer towards an off-boresight target; `--quad` keeps one polarimetric quad
(default: every channel); `--window` is the frequency window of the range profiles; `--out` is the
parent output folder (default: the project's `data/processed/bistaticity`, whatever the current
directory). Every export gets its own subfolder named after the `.npz` file, created if needed, so
analysing other data sets never overwrites earlier results (re-running the same export does). Its
outputs make no assumption about the target or the scene: no reference target, no thresholds or
"ideal" lines, no classification. It writes:

* `signatures_vs_sweep.png` and `signatures_vs_sweep.csv`: the signatures of §4.5 against the
  swept variable; per-frequency quantities as a line at the centre frequency with a band spanning
  their minimum and maximum over the frequency sweep (not a confidence interval), band-averaged
  ones as a line. Phase differences are circular: their band is the line plus the min–max of each
  phase's deviation from it, wrapped to ±180°, drawn on a [−180°, 180°] axis and repeated 360° up
  and down (a phase scattered around ±180° would otherwise span the whole circle). Bounded
  quantities use their whole domain (H, A ∈ [0, 1], ᾱ ∈ [0°, 90°], Pauli fractions up to 0 dB),
  so noise is not zoomed into structure; the other axes fit the data, bands included;
* `signatures_vs_frequency.png`: every per-frequency signature as a frequency × swept-variable map,
  the unsummarized data behind the bands, which shows whether a wide band comes from a narrow
  null, a slope or ripple across the band, or the swept variable. Each panel has its own colour
  scale over its data (ratios symmetric about 0 dB); phases are shown as their circular mean plus
  the wrapped deviation, on a cyclic scale when the deviations reach around the circle;
* `range_profiles.png`: the beamformed range profiles of the total power and of every
  polarimetric channel towards the target;
* `azimuth_elevation_maps.png`: total power and Pauli powers at the target range, only for arrays
  with several virtual channels per polarization (a cluster has no angular resolution).

Both `sbr_results` and the analysis warn when several values of the swept variable give
identical data, e.g. an SBR+ variable that does not reach the linked antenna (a stale link
re-created from scratch fixed one such case); the data are still written and analysed.

Off boresight, the steering needs the virtual channel positions; they come from the builder's port
table in the export (a hand-built design such as the horn cluster has none, so it is analysed at
boresight, where every steering vector is 1). The port table holds them in the array frame
(boresight +z, H along x, V along y), and the processor turns them into its scene frame (§4.3),
so the steering angles mean the same whatever the frame.

### 4.3 The processing chain

`MimoPolarimetryProcessor` keeps the exported channels S_c(f_n), c = S(Rx port, Tx port), for every
variation, grouped into the polarimetric channels q ∈ {HH, HV, VH, VV} (receive, transmit), each
with its virtual channels p (Tx/Rx pairs). Every step is a method; `raw()` returns the export
untouched.

| Step | Method / function | What it does |
| --- | --- | --- |
| (a) | `from_npz`, `complex_from_db_deg` | The export holds complex S directly (AEDT's real and imaginary parts, exact). Magnitude/phase data (reports, measurements) rebuild as S = 10^(M/20) e^{jπΦ/180}. Builder port signs are removed here. |
| (b) | `frequency_window` | Taper w_n over the frequency samples (Hann by default). |
| (c) | `range_profiles` | x_c[m] = IFFT(w_n S_c(f_n)) / Σw, zero-padded (4× by default), so a point target of amplitude A peaks at \|A\| at its range. Range bins are c·m/(2 N Δf), up to c/(2Δf). The phase refers to the centre frequency. |
| (d) | `spatial_weights` | Optional separable taper (`hann`, `hamming`, `blackman`) over each polarimetric channel's virtual aperture, for sidelobe control in real MIMO arrays; off by default, and pointless with one virtual channel per polarization. |
| (e) | caller | Azimuth/elevation grids, azimuth α from boresight towards H and elevation ε towards V. In the default ISO 8855 scene frame (`boresight="x"`: x forward, y left = H, z up = V), r = (cos ε cos α, cos ε sin α, sin ε); `boresight="z"` keeps the HFSS array frame (boresight +z, x = H, y = V), r = (cos ε sin α, sin ε, cos ε cos α). Positions are in the same frame: the port table's (array frame) are turned into it, while `positions=` are taken as given. Both frames give the same products. |
| (f) | `beamform` | b_q(m, α, ε) = Σ_p w_p x_q[m, p] a_p*(α, ε) / Σ_p w_p with a_p = e^{+j k₀ v_p·r} and the virtual position v_p = p_Tx + p_Rx. |
| (g) | `range_azimuth_maps`, `azimuth_elevation_maps`; `total_power`, `channel_powers`, `pauli_powers` | Maps of the mean power over the polarimetric channels, of every channel's power relative to one common reference (the total-power peak), and of the Pauli powers \|k_i\|². |

* **Sign of the steering vector.** AEDT uses e^{+jωt}, so a path L adds e^{−jkL} and a far-field
  target in direction r gives every virtual channel e^{+jk v_p·r} (times common terms). The
  matched beamformer therefore multiplies by the *conjugate* of a_p = e^{+jk₀ v_p·r}. With the
  opposite time convention the conjugate moves to the steering vector itself.
* **No Doppler.** An SBR+ frequency sweep of a static scene has no slow-time axis, so the maps are
  range–azimuth and azimuth–elevation, not range–Doppler.

`CanonicalTargetProcessor` adds `target_matrices(azimuth, elevation)`: every polarimetric channel
beamformed towards the target in the frequency domain (exact k per frequency), giving S[rx, tx]
(variation × frequency × 2 × 2). At boresight, or for a cluster with one virtual channel per
polarization, this is the plane-wave compensation of each channel towards the target. On top come
the target-agnostic `signatures` of §4.5. Interpretations are explicit opt-ins that the analysis
never uses: `fidelity(matrices, reference)` with `REFERENCE_SPHERE`, `REFERENCE_TRIHEDRAL` or
`reference_dihedral(roll)`, the `bounce_parity` classifier, and `calibrate_with_sphere`, which can
only assume equal cross-pol channel gains, because a sphere does not reveal them (see its
docstring).

```python
from pomaar.simulator.polarimetry_processor import CanonicalTargetProcessor

processor = CanonicalTargetProcessor.from_npz("horn_cluster_sphere_boresight_3m.npz")
matrices = processor.target_matrices()                         # boresight
features = processor.signatures(matrices=matrices)
range_m, profiles = processor.range_profiles()                 # raw channel range profiles
```

**Quads and port tables.** Channels are read from their port names, both hand-built
(`HornCluster1_RxH`) and builder ones (`Linked1_Port_Rx_1_H`); every channel takes part by default,
and `quad=` (a name from `available_quads(path)`, e.g. `Rx_1/Tx`, or a dict of four ports) keeps
one quad. For an array built by `hfss_array_builder`, `sbr_results` attaches the builder's port
table; the processor then removes the port field signs (`correct_signs=False` keeps them), so that
every channel refers to H and V (the array frame's +x and +y), and takes the virtual positions
from the port phase centres, turned into its scene frame. Those are build-time values, so a sweep
of a layout variable still needs `positions=` per variation. `positions` (port name or suffix →
metres, in the scene frame: ISO by default, the array frame with `boresight="z"`) also places
hand-built designs.

`sbr_results.export_sbr_parametric_sweep` stores the following in the `.npz`:

* `s` (variation × channel × frequency, complex) and `frequencies_hz`;
* `sweep_values`;
* `expressions`, i.e. the AEDT names `S(<Rx port>,<Tx port>)`;
* the nominal design variables;
* for an array built by `hfss_array_builder`, its port table (`metadata["ports"]`): per excitation
  the builder port, role, polarization, field sign and phase centre (mm). `sbr_results` finds the
  `<ArrayDesign>_ports.json` that matches the excitations next to the project, or takes
  `--port-table <file>`.

Like the other AEDT scripts, `sbr_results` attaches to the running GUI session and never starts an
AEDT.

> [!NOTE]
> If `get_solution_data` returns *"No Data Available"*, the SBR+ results were invalidated (e.g.
> by an edit of the linked source design). Re-solve the SBR+ design; the exporter then works.

**Checks.** `sbr_results` checks the output folder before it attaches to AEDT (relative paths
resolve against the current directory, which the `ansys-vnc` wrapper keeps). It refuses to write
an export in which some variations lack data at some frequencies: AEDT returns all variations on
the union of their frequency grids, with NaN where a variation was solved on another grid (e.g. a
re-solved nominal next to old parametric results) or not at all. The error lists every
incomplete variation with its number of frequencies and their step; re-solve the parametric sweep
and export again. The processor (`check_complete_sweep`) and `bistaticity_analysis` apply the
same check on load, and the analysis checks its input file and output folder first.

**Convention.** Scattering matrices are indexed `S[..., rx, tx]` with 0 = H and 1 = V, so
`S_HV` = receive H, transmit V = AEDT `S(...RxH, ...TxV)`.

### 4.4 What the exported S-parameters contain

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
at every delay. An earlier version of this study range-gated the total solution; the leakage
contaminated the co-pol channels at 3.41λ and 5.41λ and, inside the gate, the cross-pol channels:
gated total-field cross-pol came out at −30…−41 dB at 79 GHz, whereas the scattered solution
gives −36…−57 dB (§4.7). That gating chain has been removed.

> [!IMPORTANT]
> **Remove the coupling at its source.** Export `<setup> : <sweep>_Scattered` (the default).
> The processor and the analysis assume it and do no background subtraction of their own.

### 4.5 Signatures (panels of `signatures_vs_sweep.png`)

The signatures are computed from `target_matrices`, i.e. the channels beamformed towards the
steering direction (at boresight: as exported). None of them refers to a particular target.

| # | Signature | Key in `signatures()` | Per |
| --- | --- | --- | --- |
| 1 | Channel returns \|S_pq\|² | `power_<pq>_db` | frequency |
| 2 | Span ‖S‖²_F | `span_db` | frequency |
| 3 | Co-pol ratio \|S_VV/S_HH\| | `copol_ratio_db` | frequency |
| 4 | Co-pol phase difference ∠(S_VV S_HH*) | `copol_phase_difference_deg` | frequency |
| 5, 6 | Co- to cross-pol ratio per Tx pol, the two receivers of one transmission (as in a weather radar's linear depolarization ratio): magnitude \|S_co / S_cross\| and phase difference ∠(S_co S_cross*), with S_HH/S_VH for Tx H and S_VV/S_HV for Tx V | `co_minus_cross_tx<h\|v>_db`, `co_cross_phase_tx<h\|v>_deg` | frequency |
| 7, 8 | Cross-pol ratio S_HV/S_VH: magnitude and phase difference | `reciprocity_db`, `reciprocity_phase_deg` | frequency |
| 9 | Normalized asymmetry \|S_HV − S_VH\| / ‖S‖_F | `asymmetry_db` | frequency |
| 10 | Bistatic Pauli power fractions of k = [HH+VV, HH−VV, HV+VH, j(HV−VH)]/√2 | `pauli_fractions_db` | frequency |
| 11, 12 | Cloude–Pottier entropy H and anisotropy A = (λ₂−λ₃)/(λ₂+λ₃) (one panel, both in [0, 1]), and mean α of the bistatic 4×4 coherency matrix | `entropy`, `anisotropy`, `mean_alpha_deg` | variation (averaged over the band) |
| (CSV) | Pauli angle of the single matrix, cos α = \|k₁\|/‖k‖ | `alpha_deg` | frequency |

Notes on the table:

* **Bistatic 4×4 only.** The monostatic 3×3 coherency matrix symmetrizes HV and VH, i.e. assumes
  reciprocity; the analysis does not use it. k₄ vanishes for any monostatic reciprocal
  measurement, so its fraction directly shows the asymmetry of the measured matrix.
* **Cloude–Pottier over the band.** For one deterministic target, the coherency matrix is averaged
  over the frequency samples; H > 0 then means the signature changes across the band.
* **Phase reference.** All phases are relative between channels. The absolute arg S_HH is path
  length (geometry), not polarimetry.
* **Reading them (not done by the code).** An odd-bounce scatterer (sphere, trihedral, plate)
  shows S_VV ≈ S_HH (0 dB, 0°), all power in k₁ and ᾱ ≈ 0°; an even-bounce one (dihedral with its
  edge along H or V) S_VV ≈ −S_HH (0 dB, ±180°), k₂ and ᾱ ≈ 90°; a 45°-rolled dihedral is purely
  cross-polar (k₃). These are comparisons the reader makes, or opt-ins such as
  `CanonicalTargetProcessor.fidelity`.

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
helix scattering. Values computed with `plane_wave_phase` for this cluster (an earlier version of
`bistaticity_analysis` plotted them as `linear_term_study.png`):

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
for its own Tx/Rx pair and the look direction. Beamforming does exactly this, as long as the
steering vectors use each polarization channel's own virtual positions, not one shared array; this
is how `beamform` and `target_matrices` work. Better still, divide by `spherical_wave_phase` for the focal point, which also removes the
Fresnel (quadratic) term. **Polarimetric errors caused by the linear term are therefore a
processing error, not a limit of the configuration.**

**What survives compensation** sets the real bistaticity limit:

1. **Fresnel residual, when only plane-wave compensation is used** (as in the beamforming). Elements on a circle around the centroid make it common-mode at
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
4. **Coupling and multipath.** Direct coupling is removable (scattered solution or background
   subtraction).
   Interactions between the target and the antenna structure are absent from this SBR+ model
   unless the structure is included in the scene.

### 4.7 First results: boresight sphere at 3 m

> **Caveat.** These exports linked the horn cluster's interpolating sweep, so their frequency
> dependence carries interpolation ripple, with slope breaks at the sweep's basis frequencies
> (Module 2, "Only discrete sweeps are linked"). The band-wide ranges below include it; the
> 79 GHz values (a basis frequency) are less affected. Re-solve with a linked discrete sweep
> before relying on them.

From the original 201-point export (8 spacings, 50 MHz steps; since replaced in
`data/interim/bistaticity/horn_cluster_sphere_boresight_3m.npz` by a 1001-point, 20-spacing sweep).
The values are at 79 GHz, with H, A and ᾱ taken over the band; β = 2·atan(d/2R) is the co-polar
bistatic angle of this scene:

| d (λ₀) | β (°) | \|S_VV/S_HH\| (dB) | ∠(S_VV S_HH*) (°) | Co − cross, Tx H (dB) | Co − cross, Tx V (dB) | \|S_HV/S_VH\| (dB) | H (10⁻⁴) | A | ᾱ (°) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.41 | 0.11 | −0.08 | −0.4 | 41.2 | 36.5 | 4.7 | 7.7 | 0.89 | 0.74 |
| 2.41 | 0.18 | −0.07 | 0.8 | 40.6 | 57.5 | −16.9 | 3.9 | 0.46 | 0.50 |
| 3.41 | 0.26 | −0.14 | 0.6 | 40.3 | 40.5 | −0.4 | 5.4 | 0.67 | 0.40 |
| 4.41 | 0.33 | −0.19 | 0.8 | 41.1 | 48.5 | −7.6 | 4.5 | 0.84 | 1.02 |
| 5.41 | 0.41 | 0.00 | 0.6 | 42.9 | 44.2 | −1.2 | 2.3 | 0.82 | 0.52 |
| 6.41 | 0.48 | −0.12 | −0.8 | 39.1 | 47.9 | −8.9 | 3.0 | 0.79 | 0.71 |
| 7.41 | 0.56 | 0.15 | 0.1 | 43.5 | 41.4 | 2.2 | 1.6 | 0.47 | 0.56 |
| 8.41 | 0.63 | −0.02 | −0.1 | 44.2 | 42.8 | 1.4 | 2.3 | 0.48 | 0.33 |

* **Co-pol.** Within −0.28…+0.20 dB and −1.3…+2.7° over the full 74–84 GHz band, with no
  bistatic trend up to β = 0.63°, as expected for a sphere. Absolute |S_HH| is −97.6…−98.1 dB.
* **Cross-pol.** The cross-pol returns sit 33…73 dB below the co-pol ones over the band (36…57 dB
  at 79 GHz), with no trend in β.
  It varies non-monotonically from one spacing to the next, which points to the elements'
  embedded-pattern cross-pol (mutual coupling to the orthogonal neighbour changes with spacing;
  item 3 of §4.6) rather than to the bistatic geometry.
* **HV versus VH.** They differ by up to 17 dB in magnitude and by arbitrary phase. With cross-pol
  produced by different antenna pairs, nothing forces them to agree. The antisymmetric k₄
  fraction (−40…−52 dB at 79 GHz) therefore reflects antenna diversity, not bistaticity.
* **Overall.** The single-matrix Pauli angle stays ≤ 1.5° over the band and H < 10⁻³ at every
  spacing; both are set by that cross-pol. A (0.46…0.89) only says how the little non-k₁ power
  splits between the second and third eigenvalues, so it is not meaningful at this H.
* **Interim conclusion.** Up to 8.41λ at 3 m (β ≤ 0.63°), the boresight sphere shows no
  polarimetric breakdown. The configuration's polarimetric floor (cross-pol about 40 dB below
  co-pol) is set by the antennas. To find where bistaticity starts to matter, use targets with narrow bistatic
  lobes, shorter range, or off-boresight positions (§4.8).
* **Total-field results.** An earlier, now removed, gated total-field analysis flagged 3.41λ and
  5.41λ as contaminated. It agreed on co-pol at the clean spacings, but overstated the cross-pol
  (−30…−41 dB) because coupling leaked inside the gate.

### 4.8 Next steps

* **Finer frequency steps.** Halving Δf, which the SBR+ sweep copies from the array's discrete
  sweep and that in turn from the unit cell's sweep, doubles the unambiguous range of the range
  profiles and maps, at one more full-wave array solve per added frequency.
* **Cross-pol cross-check.** Check the cross-pol level against the HFSS embedded patterns:
  M_VH ≈ e_rx,V · S_ideal · e_tx,H at θ ≈ β/2 for every element.
* **Off-boresight targets.** Put the target off boresight (θ_y scan) to observe the linear term in
  simulation, and verify that per-channel compensation removes it.
* **Dihedral and trihedral.** `sbr_simulator --target dihedral|trihedral` builds both; use
  `targetRoll` for the edge parallel and perpendicular to each baseline and for the 45° dihedral
  (reciprocity), and a range R ≈ 1 m so that d·L/(λR) reaches order 1.
* **More products.** Further map quantities (e.g. coherency/decomposition maps) as functions of the
  processor; per-variation virtual positions for sweeps of a layout variable (the port table holds
  build-time values); a Doppler axis once a scene has slow time (SBR+ Range-Doppler setups or a
  moving target).