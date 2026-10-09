# Polarimetric MIMO Simulation Sub-Project

This sub-project automates a three-stage simulation pipeline linking **Ansys HFSS (full-wave)**, **Ansys HFSS SBR+ (ray-tracing)**, and **Python DSP post-processing** for polarimetric MIMO automotive radar array analysis.

---

## 1. Prerequisites & Environment Setup

The pipeline requires **Ansys Electronics Desktop (AEDT) 2025.2 or newer** and a Python environment with PyAEDT (`pyaedt`), provided by this repository's `uv` environment (`uv sync`).

### Option A: Running in the AEDT container (Linux, recommended)
On the Fedora workstation, AEDT runs in the `ansys-vnc` container (see `~/Repositories/ansys-vnc`). Python scripts must run **inside** that container: from the host, PyAEDT cannot see the container's gRPC session and would silently start a native AEDT instead (the builder refuses to do so). The container mounts `~/Repositories` and the host's `uv` Python installations, so the project's `.venv` works there unchanged:
```bash
ansys-vnc start                      # once; then open the GUI from the VNC desktop or `ansys-vnc aedt`
cd ~/Repositories/pomaar
ansys-vnc exec .venv/bin/hfss_array_builder <path_to_project.aedt> <unit_cell_design_name> [layout.yaml]
```
`ansys-vnc exec` keeps the current directory and pins the process to the compute cores. The builder connects to the running GUI session's gRPC port (found automatically; `--port` overrides it) and starts a non-graphical AEDT in the container if none is listening.

### Option B: Running natively (Windows / supported Linux)
Ensure that:
1. `ANSYSEM_ROOT252` (or your corresponding AEDT installation path) is set in your environment variables.
2. The `pomaar` package dependencies are installed (`uv sync`).
3. On Linux outside a container, set `POMAAR_ALLOW_NATIVE_AEDT=1` to let the builder start or connect to AEDT.

---

## 2. Running the Pipeline

### Module 1: HFSS Full-Wave Array Synthesis (`hfss_array_builder.py`)
This script automates the creation of a planar MIMO array on a single contiguous PCB board by copying, replicating, and boolean-cutting template geometries from an isolated unit-cell element design.

To run (prefix with `ansys-vnc exec .venv/bin/` in the container setup):
```bash
hfss_array_builder <path_to_project.aedt> <unit_cell_design_name> [layout.yaml]
```
* **Example:**
  ```bash
  hfss_array_builder ~/Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt "ApertureCoupledPatch"
  ```
* **What it does:** Creates/updates a target design named `ApertureCoupledPatchMimoArray`, syncs design variables, sizes the PCB board to perfectly fit feedlines, duplicates Tx and Rx elements according to grid positions, pre-enables post-processing, assigns wave ports (lumped-port fallback), excites the first Tx port for far-field reports, draws a $\lambda_0/4$ vacuum radiation Airbox, and sets up Theta-Phi far-field spheres in `PhaseCentreCS` with realized-gain pattern and directivity-based XPD reports. The far-field phase reference is set per port to each element's CS (`--phase-reference-mode per-port`, the default), so that a linked SBR+ antenna launches every port's rays from its own element phase centre; `single` references all ports to the array centroid `PhaseCentreCS` instead, which is only valid for targets in the far field of the whole array.

### Module 2: SBR+ Target Solver (`sbr_simulator.py`)
This script builds an SBR+ design next to a solved HFSS antenna design of the same project: it places a canonical PEC target, links the HFSS design as the antenna source, pairs every Tx port with every Rx port, copies the frequency sweep, creates quick-look S-parameter reports and optionally solves.

To run (prefix with `ansys-vnc exec .venv/bin/` in the container setup; `uv sync` installs the script):
```bash
sbr_simulator <path_to_project.aedt> <hfss_design_name> [--target dihedral] [--size 40mm] [--distance 4meter] [--map-variables VAR ...] [--field-type farfield] [--simulate | --no-simulate]
```
* **Example:**
  ```bash
  sbr_simulator ~/Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt DualPolHornCluster --target sphere --map-variables copolarSpacingLambda
  ```
* **What it does**, in the design `Linked<hfss_design_name>` (`--design` overrides; an existing one is only replaced after confirmation or `--overwrite`):
  1. Creates `AntennaCS` at the global origin and places the linked antenna in it (boresight +z).
  2. Creates `TargetCS` in `AntennaCS` at the spherical position `X = targetDistance cos(targetAzimuth) sin(targetElevation)`, `Y = targetDistance sin(targetAzimuth) sin(targetElevation)`, `Z = targetDistance cos(targetElevation)`; `targetElevation` is therefore the polar angle off boresight. The CS's Euler angles ZYZ are `(Phi, Theta, Psi) = (targetAzimuth, targetElevation, targetRoll - targetAzimuth)`: the target's local z axis points radially away from the antenna without twisting, and `targetRoll` is a pure rotation about the line of sight (e.g. a 45° dihedral). Edit the TargetCS angles directly for other orientations.
  3. Creates the target in `TargetCS`, facing the antenna along local -z, sized by `targetSize` (`--size`) and made PEC by a Perfect E boundary (sheets carry no material); its SBR+ surface roughness and height standard deviation are 0 unless `--roughness` / `--height-deviation` are given:
     * sphere: PEC solid of radius `targetSize`, centred at `TargetCS`;
     * dihedral: two square `targetSize` × `targetSize` sheets, i.e. `targetSize` is the seam length and the sheet depth. The seam lies on the local x axis through `TargetCS` (parallel to H at zero roll), which is also the phase reference of the double bounce. One sheet is drawn, tilted by 45° and copied by a 180° rotation about the symmetry axis (local z);
     * trihedral: three right-isosceles triangular sheets with edge length `targetSize`, the corner at `TargetCS` and the symmetry axis pointing at the antenna. One sheet is drawn and copied by 120° rotations about local z.

     The defaults (40 mm at 4 m) keep all three in their own far field at 79 GHz (2D²/λ ≤ 3.4 m); the script warns when the chosen range is shorter. `--target none` builds the same scene without a target; its solve is the background (the direct Tx→Rx coupling between the proxy ports), a cross-check of the incident solution below (use `--design` to keep both designs).
  4. Links the HFSS design (solution `--source-solution`, default its first frequency sweep) with *Simulate source design as needed*, *Model Visualization* of all its objects except the radiation boundary objects (`RadiatingSurface`, the builder's `Airbox`), and far-field (standard) or near-field (bumper/radome integration) results, asked interactively unless `--field-type` is given. `--map-variables` maps independent HFSS design variables to SBR+ design variables of the same name, which are created from the HFSS value when missing, so they can be swept in the SBR+ design.
  5. Selects Tx/Rx: ports whose names contain `Tx`/`Rx` and end in `H`/`V` (e.g. `TxH`, `Port_Tx_1_V`); every Tx transmits to every Rx.
  6. Creates `Setup` (4 rays per wavelength, 5 bounces, PTD + UTD, creeping waves) whose sweep `Sweep` copies the linked HFSS sweep.
  7. Creates quick-look reports. SBR+ splits the solution `Setup : Sweep` (total) into `Setup : Sweep_Incident`, the direct Tx→Rx coupling between the proxy ports, and `Setup : Sweep_Scattered`, the target return (total = incident + scattered). The coupling dominates the total at short Tx/Rx spacings, so target signatures must be read from `Setup : Sweep_Scattered`. Reports: `S-parameters dB` (total) and `S-parameters dB (scattered)` for all Tx/Rx pairs, and, per dual-polarised Tx/Rx element pair from the scattered solution, `Co-pol VV/HH` (`S_VV/S_HH` in dB and degrees; ideally 0 dB for both, with the dihedral's phase 180° away from the sphere's) and `Cross-pol isolation` (`S_HV/S_HH`, `S_VH/S_VV` in dB).
  8. Asks whether to solve (`--simulate` / `--no-simulate` decide silently). Solving also solves the HFSS design where its linked solution is missing.

---

## 3. Unit-Cell Design Rules (Strict Modeling Constraints)

To ensure the automated builder script runs successfully without manual alignment, any new unit-cell HFSS design **must** follow these strict modeling conventions:

> [!IMPORTANT]
> **Unit-Cell Modeling Rules:**
> 1. **PCB Layer Naming:**
>    * Substrates and ground planes must follow the pattern `L{Order}_{Type}`:
>      * Substrates: e.g., `L12_Substrate`, `L23_Substrate` (must end with `_Substrate`).
>      * Ground Planes: e.g., `L2_Ground` (must end with `_Ground`).
> 2. **Active Geometries & Ports:**
>    * Copper patches/feedlines must have unique layer names (e.g. `L1_Patch`, `L3_Trace`).
>    * Port excitation faces must be named `PortSheet` or `PortSheet{N}` (e.g., `PortSheet1`).
> 3. **Boolean Dummy Solids:**
>    * Slots, holes, or feed cutouts must be modeled as vacuum solids and named `f"{operation}_{target}"`:
>      * Example: `Subtract_L2_Ground` will automatically be replicated at every grid position and subtracted from the ground plane.
> 4. **PhaseCentreCS Coordinate System:**
>    * The unit cell should contain a Relative Coordinate System named exactly **`PhaseCentreCS`** at the radiation phase center of the unit cell.
>    * The script automatically retrieves `PhaseCentreCS` offsets to rotate/shift elements during replication. If it is missing (or you decline to reuse it), the builder fits it by weighted least squares to the co-polar far-field phase of the solved unit cell (solving the unit-cell setup first if needed); declining the fit falls back to a `[0,0,0]` offset.
