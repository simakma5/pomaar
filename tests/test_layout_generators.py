"""Runs every layout generator and checks its output against the builder's planning rules."""

from pathlib import Path
import re
import subprocess
import sys

import pytest
import yaml

from pomaar.simulator.array_templates import (
    LayoutError,
    PortSpec,
    UnitCellTemplate,
    array_design_name,
    layout_settings,
    plan_elements,
)

LAYOUTS = Path(__file__).resolve().parents[1] / "src" / "pomaar" / "simulator" / "layouts"
GENERATORS = sorted(LAYOUTS.glob("generate_*.py"))
V_CELL = UnitCellTemplate("Patch", (PortSpec("PortSheet_V", "V", (0.0, 0.0, 0.0)),))
DUAL_CELL = UnitCellTemplate(
    "Dual",
    (PortSpec("PortSheet_H", "H", (0.0, 0.0, 0.0)), PortSpec("PortSheet_V", "V", (0.0, 0.0, 0.0))),
)


def generate(script, tmp_path, *args):
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--unit-cell",
            "Patch",
            *args,
            "-o",
            str(tmp_path / "x.yaml"),
        ],
        check=True,
        capture_output=True,
        cwd=LAYOUTS,
    )
    return yaml.safe_load(next(tmp_path.glob("x*.yaml")).read_text())


@pytest.mark.parametrize("script", GENERATORS, ids=lambda p: p.stem)
def test_generators_write_bound_plannable_layouts(script, tmp_path):
    layout = generate(script, tmp_path, "--extra-unit-cell", "hv=Dual")
    settings = layout_settings(layout)
    assert settings.primary_unit_cell == "Patch"
    assert settings.unit_cells == {"HV": "Dual"}
    assert settings.ignored_keys == ()
    # Every element is placed by a numeric position and a parametric (x, y) expression
    for element in layout["elements"]:
        assert len(element["position"]) == 3 and len(element["position_expression"]) == 2
    plans, _ = plan_elements(
        layout["elements"], {"Patch": V_CELL, "Dual": DUAL_CELL}, "Patch", settings.unit_cells
    )
    assert plans


def test_generators_require_a_unit_cell():
    result = subprocess.run(
        [sys.executable, str(GENERATORS[0])], capture_output=True, text=True, cwd=LAYOUTS
    )
    assert result.returncode != 0 and "--unit-cell" in result.stderr


def test_layout_files_are_named_after_their_unit_cell():
    sys.path.insert(0, str(LAYOUTS))
    try:
        from _layout_common import layout_filename
    finally:
        sys.path.pop(0)

    class Args:
        unit_cell = "SeriesFedPatch"
        output = None

    name = layout_filename(Args, "contiguous_ula", "_4x4")
    assert name == Path("SeriesFedPatch_contiguous_ula_4x4.yaml")
    assert array_design_name("SeriesFedPatch", name) == "SeriesFedPatchContiguousUla4x4"


def test_shipped_layouts_are_bound_to_aperture_coupled_patch():
    for path in LAYOUTS.glob("*.yaml"):
        layout = yaml.safe_load(path.read_text())
        settings = layout_settings(layout)
        assert path.stem.startswith(f"{settings.primary_unit_cell}_")
        assert settings.ignored_keys == ()


def test_layout_variables_are_all_referenced():
    for path in LAYOUTS.glob("*.yaml"):
        metadata = yaml.safe_load(path.read_text())["metadata"]
        layout = yaml.safe_load(path.read_text())
        text = " ".join(
            [str(v) for v in metadata.get("board", {}).values()]
            + [str(v) for v in metadata.get("variables", {}).values()]
            + [str(x) for e in layout["elements"] for x in e["position_expression"]]
        )
        unused = [n for n in metadata.get("variables", {}) if not re.search(rf"\b{n}\b", text)]
        assert unused == [], f"{path.name}: unused variables {unused}"


def test_polarimetric_ula_needs_a_dual_port_cell(tmp_path):
    layout = generate(LAYOUTS / "generate_polarimetric_ula.py", tmp_path)
    with pytest.raises(LayoutError, match="dual-port"):
        plan_elements(layout["elements"], {"Patch": V_CELL}, "Patch", {})
