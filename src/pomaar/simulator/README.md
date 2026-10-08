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
`ansys-vnc exec` keeps the current directory and pins the process to the compute cores. The builder connects to the GUI session's gRPC port (`--port`, default `$ANSYS_VNC_GRPC_PORT` or 50051) and starts a non-graphical AEDT in the container if none is listening.

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
* **What it does:** Creates/updates a target design named `ApertureCoupledPatchMimoArray`, syncs design variables, sizes the PCB board to perfectly fit feedlines, duplicates Tx and Rx elements according to grid positions, pre-enables post-processing, assigns wave ports (lumped-port fallback), excites the first Tx port for far-field reports, draws a $\lambda_0/4$ vacuum radiation Airbox, and sets up Theta-Phi far-field spheres in `PhaseCentreCS` with realized-gain pattern and directivity-based XPD reports.

### Module 2: SBR+ Target Solver (`sbr_simulator.py`)
This script links the synthesized full-wave design as a composite antenna source in SBR+, imports target geometries (e.g., calibration spheres or complex vehicle CAD), and coordinates bistatic solves.

To run:
```bash
ansys-vnc exec .venv/bin/python -m pomaar.simulator.sbr_simulator <path_to_project.aedt> <synthesized_mimo_design_name>
```

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
