# Polarimetric MIMO Simulation Sub-Project

This sub-project automates a three-stage simulation pipeline linking **Ansys HFSS (full-wave)**, **Ansys HFSS SBR+ (ray-tracing)**, and **Python DSP post-processing** for polarimetric MIMO automotive radar array analysis.

---

## 1. Prerequisites & Environment Setup

The pipeline requires **Ansys Electronics Desktop (AEDT) 2025.2 or newer** and a Python environment with PyAEDT (`ansys-aedt-core`).

### Option A: Running in the Podman Container (Recommended)
If you are running in the `ansys_vnc_desktop` container, the environment is already pre-configured with Python and AEDT. Simply execute commands using:
```bash
podman exec -it -w /home/martin/Repositories/pomaar ansys_vnc_desktop <command>
```

### Option B: Running Natively (Windows / Linux)
Ensure that:
1. `ANSYSEM_ROOT252` (or your corresponding AEDT installation path) is set in your environment variables.
2. The `pomaar` package dependencies are installed:
   ```bash
   pip install numpy scipy matplotlib ansys-aedt-core>=1.1.0
   ```

---

## 2. Running the Pipeline

### Module 1: HFSS Full-Wave Array Synthesis (`hfss_array_builder.py`)
This script automates the creation of a planar MIMO array on a single contiguous PCB board by copying, replicating, and boolean-cutting template geometries from an isolated unit-cell element design.

To run:
```bash
hfss_array_builder <path_to_project.aedt> <unit_cell_design_name>
```
* **Example:**
  ```bash
  hfss_array_builder ~/Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt "ApertureCoupledPatch"
  ```
* **What it does:** Creates/updates a target design named `ApertureCoupledPatchMimoArray`, syncs design variables, sizes the PCB board to perfectly fit feedlines, duplicates Tx and Rx elements according to grid positions, pre-enables post-processing, assigns $50\,\Omega$ lumped ports, and draws a $\lambda_0/4$ vacuum radiation Airbox.

### Module 2: SBR+ Target Solver (`sbr_simulator.py`)
This script links the synthesized full-wave design as a composite antenna source in SBR+, imports target geometries (e.g., calibration spheres or complex vehicle CAD), and coordinates bistatic solves.

To run:
```bash
python3 -m pomaar.simulator.sbr_simulator <path_to_project.aedt> <synthesized_mimo_design_name>
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
>    * You **must** create a Relative Coordinate System named exactly **`PhaseCentreCS`** at the radiation phase center of the unit cell.
>    * The script automatically retrieves `PhaseCentreCS` offsets to rotate/shift elements during replication. If missing, it displays a warning and defaults to a `[0,0,0]` offset.
