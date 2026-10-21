# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

POMAAR (Polarimetric MIMO Arrays for Automotive Radars) is a PhD research repository built on the
cookiecutter-data-science layout. The active code is the HFSS / SBR+ simulation pipeline in
`src/pomaar/simulator/`, which has its own `README.md` with usage and strict unit-cell modelling
rules. `reports/research-proposal/` is a LaTeX document (latexmk with lualatex and biber). `data/`
is gitignored except for its folder structure.

## Commands

* `uv sync` installs the package editable into `.venv`, including the `hfss_array_builder`,
  `sbr_simulator`, `sbr_results` and `bistaticity_analysis` console scripts. Python >= 3.12.
* `uvx ruff check src` and `uvx ruff format src` lint and format (line length 100, isort with
  `force-sort-within-sections`). `hfss_array_builder.py` has about 24 pre-existing E501 warnings;
  don't add new ones.
* `uv run --with pytest python -m pytest tests` runs the pure-NumPy processor, template-planning,
  setup/sweep-helper and layout-generator tests. Nothing that talks to AEDT has tests; verify it against a running AEDT (below).
* `python src/pomaar/simulator/layouts/generate_<topology>.py --unit-cell <Design> [--help]`
  writes `<Design>_<topology>.yaml` next to the script (`-o` elsewhere). The shipped YAMLs are
  bound to `ApertureCoupledPatch`; the simulator README lists the exact commands.

### Running anything that talks to AEDT

AEDT 2025 R2 runs only inside the `ansys-vnc` Podman container (repo `~/Repositories/ansys-vnc`,
launcher `ansys-vnc` on PATH). Every PyAEDT script must run **inside** that container:

```bash
ansys-vnc exec hfss_array_builder <project.aedt> <layout.yaml> [-y | --build-only] [-f] [flags]
ansys-vnc exec sbr_simulator <project.aedt> <ArrayDesign> [--target ...] [flags]
ansys-vnc exec sbr_results <project> <SbrDesign> <out.npz> [--sweep-variable VAR]   # scattered
ansys-vnc solve <project.aedt> Design:Nominal:Setup   # headless batch solve
```

* **The project `.venv` works inside the container** because the host's uv Python is mounted
  there. `ansys-vnc exec`/`shell` put the nearest `.venv` first on `PATH` (`--no-venv` opts out).
* **Wrappers.** `ansys-vnc link --host bistaticity_analysis ~/Repositories/pomaar` puts wrappers for
  the console scripts into `~/.local/bin`: `hfss_array_builder`, `sbr_simulator` and `sbr_results`
  then run in the container from any mounted directory, and `bistaticity_analysis` runs on the host.
* **Leftover processes.** Stopping a command on the host does not stop it in the container;
  `ansys-vnc ps` lists AEDT sessions (with kind and gRPC port), solvers and Python processes, and
  `ansys-vnc kill <pid | :port>` ends them (it refuses the GUI without `--force`).
* **Never connect from the host.** PyAEDT on the host cannot see the container's gRPC session and
  silently starts a native AEDT, which crashes. `connect_live_session` (in the builder, used by all
  three AEDT scripts) raises `NativeAedtRefusedError` there unless `POMAAR_ALLOW_NATIVE_AEDT=1`.
* **All AEDT scripts attach to the running GUI session only.** `connect_live_session` finds the
  GUI's gRPC port (usually 50051); there is no `--port`, and nothing starts an AEDT
  (`NoLiveSessionError` without a GUI).
* **Don't edit the user's open project unasked.** Test against a copy, or, when the user agrees,
  build sandbox designs in it from layouts with a distinct file stem (e.g. `sandbox_*.yaml` in
  `data/interim/sandbox/`), which gives distinct array design names. The builder saves the project.

## Architecture: the simulation pipeline

There are three stages, each a class re-exported by `pomaar.simulator`:

1. **`MimoHfssBuilder`** (`hfss_array_builder.py`). It builds a full-wave MIMO array design from a
   solved unit-cell HFSS design plus a layout file. It is the bulk of the code.
2. **`SbrSimulationManager`** (`sbr_simulator.py`). It builds an SBR+ design with a canonical
   target (sphere, dihedral, trihedral) on a sphere around `AntennaCS`, links the array design as
   the antenna source (only a discrete sweep that saves fields; see "Linked far fields" below),
   sets Tx/Rx pairs, copies that sweep's ranges, creates reports and optionally solves.
   SBR+ total = `_Incident` (direct Tx/Rx coupling) + `_Scattered` (target); read targets from
   `Setup : Sweep_Scattered`.
3. **`MimoPolarimetryProcessor`** (`polarimetry_processor.py`). It is pure NumPy, needs no AEDT,
   and works from `.npz` sweeps exported by `sbr_results.py` (scattered solution by default). It
   runs a plain MIMO radar chain, one method each: frequency window, IFFT to range profiles
   (phase referred to the centre frequency), optional spatial window, beamforming with steering
   vectors a_p = e^{+jk₀ v_p·r} (applied conjugated; v_p = p_Tx + p_Rx; AEDT's e^{+jωt}), and
   range–azimuth / azimuth–elevation maps of total power, channel powers and Pauli powers.
   `raw()` forwards the export. No gating and no Doppler (static frequency sweeps).
   `CanonicalTargetProcessor` adds `target_matrices` (frequency-domain beamforming towards the
   target) and target-agnostic `signatures` (channel powers, span, co-pol ratio/phase difference,
   co−cross dB, cross-pol ratio, asymmetry, bistatic Pauli fractions, Cloude–Pottier H/A/ᾱ of
   the 4×4 band-averaged coherency matrix); `fidelity`, `bounce_parity` and the reference
   matrices are opt-ins the analysis never uses. Exports are always the scattered solution:
   neither the processor nor the analysis handles total fields or background subtraction.
   `bistaticity_analysis.py` only combines and plots those products against the swept variable.

### Layouts

`layouts/*.yaml` hold `metadata` and `elements`:

* `metadata`: `primary_unit_cell` (required; the layout is bound to its designs), optional
  `unit_cells` (`{h: PatchH}`, `{hv: DualPatch}`), `drop_unused_ports`, `phase_reference_mode`,
  design `variables` (each referenced; a test checks it) and `board` size formulas
  (`array_templates.layout_settings` reads and checks them). No frequencies: the builder takes
  the centre frequency and the frequency sampling from the primary unit cell's setup. Older layouts'
  `centre_frequency_ghz`, `bandwidth_ghz` and `topology` are reported as ignored;
* `elements`: each with role Tx/Rx, polarization h/v/hv, numeric `position` (x y z, for the
  builder's checks), a parametric `position_expression` [x, y] in design variables (for the
  placement in AEDT), `yaw` and an optional `unit_cell` design.

The layout is the builder's only description of the array (`MimoHfssBuilder(project, layout)`);
CLI flags only answer prompts. The array design is named
`<primary_unit_cell><LayoutFileStemInCamelCase>`, without a leading `<primary_unit_cell>_` in the
stem (`array_design_name`). The generators share `layouts/_layout_common.py` (`--unit-cell`,
`--extra-unit-cell`, `--drop-unused-ports`, `--phase-reference-mode`, `-o`; default file
`<UnitCell>_<topology>.yaml`). The frequency sweep is never in a layout: the builder copies the
primary unit cell's setup for its meshing (or fails), replaces its sweeps by a discrete sweep over
the unit cell's sweep frequencies, and `sbr_simulator` copies that sweep.

### Polarization templates

Polarization is a property of the unit-cell template, not of placement. `array_templates.py` (pure,
AEDT-free, tested in `tests/test_array_templates.py`) holds the rules:

* **Tags.** Port sheets must be tagged `PortSheet[<n>]_<H|V>`; an untagged sheet fails the
  polarization check. A cell has one tagged port or two (dual-port).
* **Resolution.** An element uses its `unit_cell`, else the template mapped to all its
  polarizations, else the source design turned by its yaw.
* **Check.** `plan_elements` checks that R(yaw)·native lies along x (H) or y (V) before Step 2,
  so an abort never clears a design.
* **Unused ports** of a dual-port cell become matched `Port_<label>_<P>_Load` ports, which
  `sbr_simulator` skips.
* **Naming.** Ports are `Port_<label>_<P>` (`_parse_port` pairs H/V by the shared prefix), with
  one CS per active port, `CS_` + the port name without `Port_`.
* **Signs.** Yaw 180 flips E. The builder prints the per-port signs, warns when a role's pair
  sign s_H·s_V = −1 (cross-pol negated relative to co-pol, invisible to sphere calibration), and
  writes `<TargetDesign>_ports.json` next to the project.
* **`unit_cells` keys** are polarization sets: `h`/`v` serve single-pol elements, `hv` only
  dual-polarized (`polarization: hv`) ones.

### How the stages share port information

`array_templates.py` also holds the channel helpers that all stages use: `parse_port` (the last
`Tx|Rx` counts, so instance prefixes such as `Linked1_Port_Rx_1_H` parse), `s_parameter_ports`,
`polarimetric_quads` (an H/V pair per role, either sharing an element key or the role's single H
and V ports) and `attach_port_table`.

* `sbr_simulator` skips `_Load` ports, builds its reports from the scattered solution (one
  total-field report for reference), per quad for VV/HH and cross-pol, and sign-corrects the VV/HH
  phase from the source design's `_ports.json`.
* `sbr_results` exports `Setup : Sweep_Scattered` by default and attaches the matching
  `_ports.json` (found next to the project, or `--port-table`) to the `.npz` metadata as `ports`.
* `MimoPolarimetryProcessor.from_npz` keeps every channel (or one `quad=`), removes the port
  signs (`correct_signs`) and takes the virtual positions from the port phase centres
  (`positions=` overrides, e.g. for hand-built designs). Positions are build-time values, so they
  don't follow a swept layout variable; without any, steering is exact only at boresight.

### Builder flow

`synthesize_array_in_hfss`, then `create_post_processing_reports`, then `run_solve`. The numbered
`STEP` comments in the code are the map:

* **Steps 1–3:** classify the objects of every template by name convention and plan the
  elements, create or reuse the target design (variables merged from all templates, setup from
  the source design), then size and draw the board. All templates must share the board stack,
  which sits at the source design's phase-centre height `z_ref`.
* **Step 4:** each active port gets its own relative CS at its phase centre with the element yaw.
  Each template's geometry is copied once and tagged `<obj>__<template>` (so equal names in two
  unit cells don't clash), then replicated with `duplicate_along_line`; `map_duplicates` maps the
  copies back exactly. `.clone()` (copy/paste) hangs on a headless X11/VNC display.
* **Step 5:** `<Operation>_<Target>` dummy solids are applied as booleans.
* **Step 6:** wave ports (and `_Load` ports) are assigned on `PortSheet*`, falling back to lumped
  ports.
* **Step 6b:** `PhaseCentreCS` is placed at the array centroid as the CS of the far-field spheres.
  The far-field phase reference is per port (each port CS), so SBR+ launches every port from
  its own phase centre. Only the CS origin appears to matter, not its yaw (no co/cross swap seen).
* **Step 7:** a laterally flush radiation airbox is added (a deliberate choice: wave ports must sit
  on the boundary).
* **Step 7b:** `verify_placement` checks every element against R(yaw)(p − o) + pos (its template's
  offset o, z shifted by `z_ref`) and aborts on mismatch. Step 7c reports the port polarizations.
* **Setup (`configure_simulation_setup`, after Step 7c):** `Setup` is the unit cell's setup
  (native copy-paste, else recreated) for its meshing only. Its copied sweeps are replaced by one
  discrete `Sweep` over the ranges of the unit cell's first sweep (read with `saved_setups`) that
  saves the radiated fields; a unit cell without a sweep fails. It warns when the adaptive
  frequencies stay below the top of the sweep. A discrete sweep is one full-wave solve per
  frequency.
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
* **Linked far fields are interpolated between solved frequencies.** SBR+ gets a linked antenna's
  far fields only where the HFSS design saved fields: an adaptive solution (`LastAdaptive`) has
  them at the adaptive frequencies, an interpolating sweep at its basis frequencies (bisections
  of the band). In between, SBR+ interpolates, which put ripple with slope breaks exactly at those
  frequencies into every result (seen at 77/79/81 GHz and at 75.25, 76.5, …, 82.75 GHz). Hence the
  builder's discrete `Setup : Sweep` (radiated fields saved) and `sbr_simulator`'s refusal of
  anything else (`frequency_sweeps.linked_solution_problem`). The SBR+ frequencies must be among
  the saved ones; `sbr_simulator` copies the linked sweep's ranges to ensure it.
* **PyAEDT sweep objects misreport the type.** `SweepHFSS.__init__` (PyAEDT 1.1) sets every
  loaded sweep's `props["Type"]` to `Interpolating`, and setups and sweeps are read from the saved
  `.aedt`, which may predate GUI edits. Read them with the builder's `saved_setups(app)` (saves
  the project, then parses the file with `load_entire_aedt_file`). `oanalysis.GetSweeps` and the
  object API (`odesign.GetChildObject("Analysis")`, `GetPropValue("Type")`) are live but lack the
  save-fields flags.
* **Mapped variables of a linked antenna.** `InsertNativeComponent` accepts a mapping
  (`Params` `name='name'`, instance `GeometryParameters`) whose name is no parameter of the
  component definition. Such a link keys its per-instance caches (`InstanceBoundingBoxMap`) by no
  parameter: the variable resizes only the `_Footprint` box, not the visualized antenna, and a
  horn link made that way gave identical SBR+ results at every swept value. `link_antenna` therefore also declares each mapped
  variable in `GeometryDefinitionParameters` as `"VariableProp:=", [name, "D", "", value]`, as the
  GUI does. A correct link's saved definition shows `LocalVariables` with that `VariableProp` and
  one bounding-box `Key` per value used.
* **Per-variation solutions.** `existing_analysis_sweeps` also lists `Setup - <variation> : …`
  entries (e.g. from SBR+ "simulate source design as needed"); only names whose setup is in
  `setup_names` are the design's own.

## Conventions

* **Analysis outputs are unbiased.** Figures and CSVs show quantities that assume nothing about
  the target or the scene: no reference target or fidelity to one, no "ideal" values, threshold
  lines or classifications, no titles naming the target. Interpretation belongs in the README
  text or in explicit opt-in functions (e.g. `CanonicalTargetProcessor.fidelity`). Bands in the
  signature plots are the min–max over the frequency sweep, not confidence intervals.

* British spelling with Oxford -ize in code and docs (centre, metre, polarization). Keep AEDT,
  PyAEDT and matplotlib names as they are (`SetPhaseCenterPerPort`, `analyze()`, the `meter`
  unit token, `color=`). CLI flags accept American aliases, and names read from AEDT designs are
  looked up in both spellings (`_spelling_variants`/`_find_by_spelling` in the builder).

* Commit messages: `simulator: <past-tense summary>`, e.g. `simulator: fixed far-field
  definitions, ...`, or plain past tense for repo-wide changes.
* `*_notes.md` files (e.g. `src/pomaar/simulator/implementation_notes.md`) are gitignored working
  notes.
* The top-level `README.md` tree lists `io.py`, which does not exist. The core package is
  `config.py` (paths from `PROJ_ROOT`), `array_synthesizer.py` (MIMO topology synthesis/plots) and
  stub `dsp.py`/`plots.py`.
