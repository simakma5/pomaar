# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

POMAAR (Polarimetric MIMO Arrays for Automotive Radars) is a PhD research repository built on the
cookiecutter-data-science layout. The active code is the HFSS / SBR+ simulation pipeline in
`src/pomaar/simulator/`, which has its own `README.md` with usage and strict unit-cell modelling
rules. `reports/research-proposal/` is a LaTeX document (latexmk with lualatex and biber). `data/`
is gitignored except for its folder structure.

## Commands

* `uv sync` installs the package editable into `.venv`, including the `hfss_array_builder` and
  `sbr_simulator` console scripts. Python >= 3.12.
* `uvx ruff check src` and `uvx ruff format src` lint and format (line length 100, isort with
  `force-sort-within-sections`). `hfss_array_builder.py` has about 39 pre-existing E501 warnings;
  don't add new ones.
* `uv run --with pytest python -m pytest tests` runs the pure-NumPy processor tests. Nothing that
  talks to AEDT has tests; verify it against a running AEDT (below).
* `python src/pomaar/simulator/layouts/generate_<topology>.py [--help]` regenerates a layout YAML
  next to the script.

### Running anything that talks to AEDT

AEDT 2025 R2 runs only inside the `ansys-vnc` Podman container (repo `~/Repositories/ansys-vnc`,
launcher `ansys-vnc` on PATH). Every PyAEDT script must run **inside** that container:

```bash
ansys-vnc exec hfss_array_builder <project.aedt> <UnitCellDesign> [layout.yaml] [flags]
ansys-vnc exec sbr_simulator <project.aedt> <ArrayDesign> [--target ...] [flags]
ansys-vnc exec python -m pomaar.simulator.sbr_results <project> <SbrDesign> <out.npz> \
    --setup "Setup : Sweep_Scattered" [--sweep-variable VAR]   # export for the processor
ansys-vnc solve <project.aedt> Design:Nominal:Setup   # headless batch solve
```

* **The project `.venv` works inside the container** because the host's uv Python is mounted
  there. `ansys-vnc exec`/`shell` put the nearest `.venv` first on `PATH` (`--no-venv` opts out).
* **Leftover processes.** Stopping a command on the host does not stop it in the container;
  `ansys-vnc ps` lists AEDT sessions (with kind and gRPC port), solvers and Python processes, and
  `ansys-vnc kill <pid | :port>` ends them (it refuses the GUI without `--force`).
* **Never connect from the host.** PyAEDT on the host cannot see the container's gRPC session and
  silently starts a native AEDT, which crashes. The builder raises `NativeAedtRefusedError` there
  unless `POMAAR_ALLOW_NATIVE_AEDT=1`.
* **The builder attaches to the running GUI session.** It discovers the GUI's gRPC port (usually
  50051) and starts a non-graphical AEDT only if none is running.
* **Don't edit the user's open project.** Test against a copy, e.g. a sandbox `.aedt`, or a
  separate container/session.

## Architecture: the simulation pipeline

There are three stages, each a class re-exported by `pomaar.simulator`:

1. **`MimoHfssBuilder`** (`hfss_array_builder.py`). It builds a full-wave MIMO array design from a
   solved unit-cell HFSS design plus a layout file. It is the bulk of the code.
2. **`SbrSimulationManager`** (`sbr_simulator.py`). It builds an SBR+ design with a canonical
   target (sphere, dihedral, trihedral) on a sphere around `AntennaCS`, links the array design as
   the antenna source, sets Tx/Rx pairs, copies its sweep, creates reports and optionally solves.
   SBR+ total = `_Incident` (direct Tx/Rx coupling) + `_Scattered` (target); read targets from
   `Setup : Sweep_Scattered`.
3. **`MimoPolarimetryProcessor`** (`polarimetry_processor.py`). It is pure NumPy, needs no AEDT,
   and works from `.npz` sweeps exported by `sbr_results.py` (export `Setup : Sweep_Scattered`).
   Its module-level functions handle range gating, geometric phase compensation and polarimetric
   signatures (bistatic Pauli, fidelity, Cloude α/entropy). The class chains them, adding sphere
   calibration and noise injection. `bistaticity_analysis.py` plots the signatures against Tx–Rx
   spacing; it is still specific to the boresight-sphere horn-cluster study (README §4).

### Layouts

`layouts/*.yaml` hold `metadata` and `elements`:

* `metadata`: centre frequency, bandwidth, design `variables` and `board` size formulas;
* `elements`: each with role Tx/Rx, polarization h/v, numeric `position`, a parametric
  `position_expression` in design variables, and `yaw`.

The target design is named `<UnitCell><LayoutFileStemInCamelCase>`, or `<UnitCell>MimoArray`
without a layout.

### Builder flow

`synthesize_array_in_hfss`, then `create_post_processing_reports`, then `run_solve`. The numbered
`STEP` comments in the code are the map:

* **Steps 1–3:** classify unit-cell objects by name convention, create or reuse the target design
  (variables and setup inherited from the unit cell), then size and draw the board.
* **Step 4:** each element gets its own relative CS at its position with its yaw. Template geometry
  is copied once, then replicated with `duplicate_along_line`. `.clone()` (copy/paste) hangs on a
  headless X11/VNC display.
* **Step 5:** `<Operation>_<Target>` dummy solids are applied as booleans.
* **Step 6:** wave ports are assigned on `PortSheet*`, falling back to lumped ports.
* **Step 6b:** `PhaseCentreCS` is placed at the array centroid as the CS of the far-field spheres.
  The far-field phase reference is per port (each element CS), so SBR+ launches every port from
  its own phase centre. Only the CS origin appears to matter, not its yaw (no co/cross swap seen).
* **Step 7:** a laterally flush radiation airbox is added (a deliberate choice: wave ports must sit
  on the boundary).
* **Step 7b:** `verify_placement` checks every element against R(yaw)(p − o) + pos and aborts on
  mismatch.
* **`PhaseCentreCS` in the unit cell** is either reused or fitted by `fit_unit_cell_phase_centre`:
  a |E|-weighted least-squares fit of the co-polar Ludwig-3 phase over θ ≤ 40°. Lateral components
  below 10 µm are snapped to 0.

### Far-field conventions (settled; don't revert)

* **Spheres.** All builder spheres are Theta-Phi in `PhaseCentreCS` with an explicit `Boresight="Z
  Axis"` and Linear polarization, so CoPolar/CrossPolar are Ludwig-3 X/Y. A non-Z boresight swaps
  co and cross.
* **Cuts.** The principal cuts use signed θ ∈ [−180°, 180°] at φ = 0° / 90°.
* **Metrics.** Patterns use realized gain. XPD uses directivity (`dB(DirCoPolar) -
  dB(DirCrossPolar)`), since XPD is identical for Dir/Gain/RealizedGain.
* **Excitation.** Edit Sources excites only the first Tx port (1 W).

### PyAEDT / AEDT pitfalls hit in this code

* **Working coordinate system.** `create_coordinate_system` makes the new CS the working CS, and
  subsequent move/rotate then act in it. Reset with `set_working_coordinate_system("Global")`
  right after creating one.
* **Inserting designs.** PyAEDT's `insert_design` re-points the app object it is called on to the
  new design, so later calls on that object silently act on the wrong design. Use the native
  `oproject.InsertDesign` / `DeleteDesign` when that app must stay where it is.
* **Re-creating a design.** Within one session, AEDT refuses to re-create a design under a
  just-deleted name until the project is saved: `InsertDesign` silently returns `None`. Rename the
  old design before deleting it, as `sbr_simulator.py` does.
* **`post.delete_report(None)`** deletes **all** reports; always pass a name.
* **`edit_sources`** sets every unlisted port to 0 W.
* **Reading far-field sweeps.** Use `get_solution_data(...)` and `full_matrix_real_imag`; columns
  follow `sd.intrinsics` key order. `get_expression_data` returns only one secondary slice.
* **Object caching.** PyAEDT caches modeler objects within a process, so re-read in a fresh process
  when verifying edits.
* **Compound expressions.** `evaluate_expression` can't evaluate them.
* **Batch-solve argument order.** `ansysedt -batchsolve` needs `-logfile` *before* `-batchsolve`,
  otherwise the `Design:Nominal:Setup` target is ignored and the whole project is solved.
* **Package imports are heavy.** `import pomaar` pulls in matplotlib, and `pomaar.simulator` pulls
  in PyAEDT. AEDT's bundled CPython (NumPy 1.x) can't import the package; use the project
  `.venv`.
* **SBR+ Tx/Rx pairs.** `SetSBRTxRxSettings` needs each Tx's Rx antennas as one comma-separated
  string; a Python list leaves the Tx/Rx matrix empty.
* **PEC sheets.** Sheets keep no material through PyAEDT; make them PEC with `assign_perfect_e`
  (in SBR+ it also takes the surface roughness and height deviation).
* **`create_rectangle` orientation.** The string `"Z"` does not mean "normal to z" (it drew in the
  XZ plane); pass `Plane.XY`.
* **SBR+ scattered solution.** PyAEDT's `create_report` rejects `Setup : Sweep_Scattered` before the
  first solve; the native `ReportSetup.CreateReport` accepts it.

## Conventions

* Commit messages: `simulator: <past-tense summary>`, e.g. `simulator: fixed far-field
  definitions, ...`, or plain past tense for repo-wide changes.
* `*_notes.md` files (e.g. `src/pomaar/simulator/implementation_notes.md`) are gitignored working
  notes.
* The top-level `README.md` tree lists `io.py`, which does not exist. The core package is
  `config.py` (paths from `PROJ_ROOT`), `array_synthesizer.py` (MIMO topology synthesis/plots) and
  stub `dsp.py`/`plots.py`.
