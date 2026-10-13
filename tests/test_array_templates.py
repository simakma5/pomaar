"""Tests of the pure planning helpers of the array builder (no AEDT needed)."""

import itertools

import pytest

from pomaar.simulator.array_templates import (
    LayoutError,
    PolarizationMismatchError,
    PortSpec,
    UnitCellTemplate,
    array_design_name,
    attach_port_table,
    built_polarization,
    cs_name_for_port,
    element_label,
    layout_settings,
    map_duplicates,
    merge_unit_cells,
    offset_expression,
    pair_signs,
    parse_polarizations,
    parse_port,
    place,
    place_xy,
    plan_elements,
    polarimetric_quads,
    port_name,
    port_sheet_polarization,
    port_table,
    s_parameter_ports,
    select_quad,
    sidecar_dict,
    sign_warnings,
)


def single_pol(design, pol, offset=(0.0, 0.0, 0.0)):
    sheet = f"PortSheet_{pol}"
    return UnitCellTemplate(
        design=design, ports=(PortSpec(sheet, pol, offset),), offset=offset, port_sheets=[sheet]
    )


def dual_port(design, offset=(0.0, 0.0, 0.0)):
    ports = (PortSpec("PortSheet_H", "H", offset), PortSpec("PortSheet_V", "V", offset))
    return UnitCellTemplate(design=design, ports=ports, offset=offset)


def element(label, role, pol, yaw, position=(0.0, 0.0, 0.0), **extra):
    return {
        "label": label,
        "role": role,
        "polarization": pol,
        "yaw": yaw,
        "position": list(position),
        **extra,
    }


# Element list of the shipped polarimetric_cluster.yaml (V-native cell, H by rotation)
CLUSTER = [
    element("Tx", "Tx", "h", 90.0, (1.33, 0.0, 0.0)),
    element("Rx", "Rx", "h", 270.0, (-1.33, 0.0, 0.0)),
    element("Tx", "Tx", "v", 180.0, (0.0, 1.33, 0.0)),
    element("Rx", "Rx", "v", 0.0, (0.0, -1.33, 0.0)),
]


def test_port_sheet_polarization():
    assert port_sheet_polarization("PortSheet_V") == "V"
    assert port_sheet_polarization("PortSheet1_h") == "H"
    assert port_sheet_polarization("PortSheet") is None
    assert port_sheet_polarization("PortSheet1") is None


def test_parse_polarizations():
    assert parse_polarizations("v") == ("V",)
    assert parse_polarizations("VH") == ("H", "V")
    assert parse_polarizations(["v", "h"]) == ("H", "V")
    assert parse_polarizations(None) == ()
    assert parse_polarizations("") == ()
    with pytest.raises(LayoutError):
        parse_polarizations("x")


def test_unit_cell_mappings():
    assert merge_unit_cells({}, {"h": "B", "v": "A"}) == {"V": "A", "H": "B"}
    assert merge_unit_cells({"V": "A"}, {"hv": "B"}) == {"V": "A", "HV": "B"}
    with pytest.raises(LayoutError):
        merge_unit_cells({"HV": "A"}, {"vh": "B"})
    with pytest.raises(LayoutError):
        merge_unit_cells({}, {"x": "B"})


def test_layout_settings():
    layout = {
        "metadata": {
            "primary_unit_cell": "Patch",
            "unit_cells": {"hv": "Dual"},
            "centre_frequency_ghz": 79,
            "topology": "cluster",
            "phase_reference_mode": "per-port",
        },
        "elements": [],
    }
    settings = layout_settings(layout)
    assert settings.primary_unit_cell == "Patch"
    assert settings.unit_cells == {"HV": "Dual"}
    assert settings.phase_reference_mode == "per_port"
    # The frequency setup comes from the unit cell; older layouts' keys are reported as ignored
    assert settings.ignored_keys == ("centre_frequency_ghz", "topology")
    assert not settings.drop_unused_ports
    with pytest.raises(LayoutError, match="primary_unit_cell"):
        layout_settings({"metadata": {}, "elements": []})
    with pytest.raises(LayoutError, match="phase_reference_mode"):
        layout_settings(
            {"metadata": {"primary_unit_cell": "P", "phase_reference_mode": "x"}, "elements": []}
        )
    assert array_design_name("Patch", "layouts/polarimetric_cluster.yaml") == (
        "PatchPolarimetricCluster"
    )
    assert array_design_name("Patch", "sandbox-contiguous_ula_4x4.yaml") == (
        "PatchSandboxContiguousUla4x4"
    )
    # A layout file named after its unit cell does not repeat it in the design name
    assert array_design_name("Patch", "layouts/Patch_polarimetric_cluster.yaml") == (
        "PatchPolarimetricCluster"
    )


def test_hv_unit_cell_serves_only_dual_polarized_elements():
    templates = {"PatchV": single_pol("PatchV", "V"), "Dual": dual_port("Dual")}
    plans, _ = plan_elements(
        [element("Rx_1", "Rx", "hv", 0.0), element("Tx_1", "Tx", "v", 180.0, (5.0, 0.0, 0.0))],
        templates,
        "PatchV",
        {"HV": "Dual"},
    )
    assert [plan.template.design for plan in plans] == ["Dual", "PatchV"]


def test_placement_rule():
    # The offset point lands on the position, whatever the yaw
    for yaw in (0.0, 37.0, 90.0, 180.0, 270.0):
        assert place_xy((0.3, -0.2), (0.3, -0.2), yaw, (5.0, 6.0)) == pytest.approx((5.0, 6.0))
    assert place_xy((1.0, 0.0), (0.0, 0.0), 90.0, (5.0, 6.0)) == pytest.approx((5.0, 7.0))
    # Common stack shift: z -> z - z_ref + pos_z
    assert place((0.0, 0.0, 0.5), (0.0, 0.0, 0.2), 0.2, 0.0, (0.0, 0.0, 1.0)) == pytest.approx(
        (0.0, 0.0, 1.3)
    )


def test_built_polarization():
    assert built_polarization("V", 0.0)[:2] == ("V", 1)
    assert built_polarization("V", 180.0)[:2] == ("V", -1)
    assert built_polarization("V", 90.0)[:2] == ("H", -1)
    assert built_polarization("V", 270.0)[:2] == ("H", 1)
    assert built_polarization("H", 90.0)[:2] == ("V", 1)
    assert built_polarization("V", 0.5)[:2] == ("V", 1)
    assert built_polarization("V", 45.0)[:2] == (None, 0)
    assert built_polarization(None, 0.0)[:2] == (None, 0)


def test_naming_round_trip():
    assert element_label("Rx_1", ("V",)) == "Rx_1_V"
    assert element_label("Rx_1_V", ("V",)) == "Rx_1_V"
    assert element_label("Rx_1", ("H", "V")) == "Rx_1"
    assert port_name("Rx_1_V", "V") == "Port_Rx_1_V"
    assert port_name("Rx_1", "H") == "Port_Rx_1_H"
    assert port_name("Rx_1", "H", index=2) == "Port_Rx_1_F2_H"
    assert port_name("Rx_1", "H", load=True) == "Port_Rx_1_H_Load"
    assert cs_name_for_port("Port_Rx_1_V") == "CS_Rx_1_V"
    assert parse_port("Port_Rx_1_V") == ("Rx", "V", "Port_Rx_1_")
    assert parse_port("Port_Rx_1_H")[2] == parse_port("Port_Rx_1_V")[2]
    assert parse_port("Port_Rx_1_H_Load") is None


def test_legacy_cluster_plan_is_unchanged():
    templates = {"PatchV": single_pol("PatchV", "V")}
    plans, warnings = plan_elements(CLUSTER, templates, "PatchV", {})
    assert warnings == []
    assert [plan.label for plan in plans] == ["Tx_H", "Rx_H", "Tx_V", "Rx_V"]
    ports = [plan.ports[0] for plan in plans]
    assert [p.name for p in ports] == ["Port_Tx_H", "Port_Rx_H", "Port_Tx_V", "Port_Rx_V"]
    assert [p.cs_name for p in ports] == ["CS_Tx_H", "CS_Rx_H", "CS_Tx_V", "CS_Rx_V"]
    assert [p.sheet for p in ports] == [f"PortSheet_V_{p.name[5:]}" for p in ports]
    assert pair_signs(port_table(plans)) == {"Tx": 1, "Rx": 1}
    assert sign_warnings(port_table(plans)) == []


def test_ula_yaws_flip_the_cross_pol_channels():
    # The yaws of the original polarimetric_ula.yaml: Rx V/H 0/90, Tx V/H 180/270
    layout = [
        element("Rx_1", "Rx", "v", 0.0, (0.0, -1.0, 0.0)),
        element("Rx_2", "Rx", "h", 90.0, (2.0, -1.0, 0.0)),
        element("Tx_1", "Tx", "v", 180.0, (0.0, 1.0, 0.0)),
        element("Tx_2", "Tx", "h", 270.0, (4.0, 1.0, 0.0)),
    ]
    plans, _ = plan_elements(layout, {"PatchV": single_pol("PatchV", "V")}, "PatchV", {})
    rows = port_table(plans)
    assert pair_signs(rows) == {"Tx": -1, "Rx": -1}
    warnings = sign_warnings(rows)
    assert len(warnings) == 2 and all("opposite signs" in w for w in warnings)


def test_pair_sign_is_not_the_vacuous_product_check():
    # s_HH s_VV == s_HV s_VH for every sign assignment, so it cannot detect anything;
    # the pair signs can.
    for s_rh, s_rv, s_th, s_tv in itertools.product((1, -1), repeat=4):
        assert (s_rh * s_th) * (s_rv * s_tv) == (s_rh * s_tv) * (s_rv * s_th)
    rows = [
        {"port": "Port_Tx_H", "role": "Tx", "polarization": "H", "sign": 1, "active": True},
        {"port": "Port_Tx_V", "role": "Tx", "polarization": "V", "sign": -1, "active": True},
        {"port": "Port_Rx_H", "role": "Rx", "polarization": "H", "sign": 1, "active": True},
        {"port": "Port_Rx_V", "role": "Rx", "polarization": "V", "sign": 1, "active": True},
    ]
    assert pair_signs(rows) == {"Tx": -1, "Rx": 1}
    assert len(sign_warnings(rows)) == 1


def test_mixed_signs_within_a_polarization_warn():
    rows = [
        {"port": "Port_Tx_1_V", "role": "Tx", "polarization": "V", "sign": 1, "active": True},
        {"port": "Port_Tx_2_V", "role": "Tx", "polarization": "V", "sign": -1, "active": True},
    ]
    assert pair_signs(rows)["Tx"] is None
    assert any("opposite signs" in w for w in sign_warnings(rows))


def test_wrong_polarization_aborts():
    templates = {"PatchV": single_pol("PatchV", "V")}
    with pytest.raises(PolarizationMismatchError, match="unit_cells"):
        plan_elements([element("Rx_1", "Rx", "h", 0.0)], templates, "PatchV", {})
    with pytest.raises(PolarizationMismatchError):
        plan_elements([element("Rx_1", "Rx", "v", 45.0)], templates, "PatchV", {})
    with pytest.raises(PolarizationMismatchError, match="dual-port"):
        plan_elements([element("Rx_1", "Rx", "hv", 0.0)], templates, "PatchV", {})


def test_untagged_template_needs_a_declared_polarization():
    legacy = UnitCellTemplate(design="Patch", ports=(PortSpec("PortSheet", None, (0, 0, 0)),))
    with pytest.raises(PolarizationMismatchError, match="tag its port sheets"):
        plan_elements([element("Rx_1", "Rx", "v", 0.0)], {"Patch": legacy}, "Patch", {})
    # Without a declared element polarization nothing is checked (legacy port names)
    plans, _ = plan_elements([element("Rx_1", "Rx", None, 0.0)], {"Patch": legacy}, "Patch", {})
    assert plans[0].ports[0].name == "Port_Rx_1"


def test_separate_templates_per_polarization():
    templates = {"PatchV": single_pol("PatchV", "V"), "PatchH": single_pol("PatchH", "H")}
    unit_cells = {"H": "PatchH"}
    layout = [
        element("Tx", "Tx", "h", 180.0, (1.0, 0.0, 0.0)),
        element("Tx", "Tx", "v", 180.0, (0.0, 1.0, 0.0)),
    ]
    plans, _ = plan_elements(layout, templates, "PatchV", unit_cells)
    assert [plan.template.design for plan in plans] == ["PatchH", "PatchV"]
    assert [plan.ports[0].name for plan in plans] == ["Port_Tx_H", "Port_Tx_V"]
    assert pair_signs(port_table(plans))["Tx"] == 1


def test_dual_port_cell():
    templates = {"PatchV": single_pol("PatchV", "V"), "Dual": dual_port("Dual")}
    unit_cells = {"H": "Dual", "V": "Dual"}
    plans, _ = plan_elements(
        [element("Rx_1", "Rx", "hv", 0.0), element("Tx_1", "Tx", "v", 0.0, (5.0, 0.0, 0.0))],
        templates,
        "PatchV",
        unit_cells,
    )
    rx, tx = plans
    assert rx.label == "Rx_1"
    assert [p.name for p in rx.ports] == ["Port_Rx_1_H", "Port_Rx_1_V"]
    assert [p.name for p in tx.ports] == ["Port_Tx_1_H_Load", "Port_Tx_1_V"]
    assert [p.active for p in tx.ports] == [False, True]
    assert tx.ports[0].cs_name is None
    dropped, _ = plan_elements(
        [element("Tx_1", "Tx", "v", 0.0)], templates, "PatchV", unit_cells, drop_unused_ports=True
    )
    assert dropped[0].ports[0].name is None


def test_explicit_unit_cell_and_duplicate_checks():
    templates = {"PatchV": single_pol("PatchV", "V"), "Other": single_pol("Other", "V")}
    plans, _ = plan_elements(
        [element("Rx_1", "Rx", "v", 0.0, unit_cell="Other")], templates, "PatchV", {}
    )
    assert plans[0].template.design == "Other"
    with pytest.raises(LayoutError, match="Unknown unit cell"):
        plan_elements([element("Rx_1", "Rx", "v", 0.0, unit_cell="Nope")], templates, "PatchV", {})
    with pytest.raises(LayoutError, match="Duplicate element labels"):
        plan_elements(
            [element("Rx_1", "Rx", "v", 0.0), element("Rx_1", "Rx", "v", 0.0, (9, 0, 0))],
            templates,
            "PatchV",
            {},
        )
    with pytest.raises(LayoutError, match="role"):
        plan_elements([element("Rx_1", "Tx", "v", 0.0)], templates, "PatchV", {})


def test_co_located_elements_warn():
    templates = {"PatchV": single_pol("PatchV", "V")}
    _, warnings = plan_elements(
        [element("Rx_1", "Rx", "v", 0.0), element("Rx_1", "Rx", "h", 270.0)],
        templates,
        "PatchV",
        {},
    )
    assert len(warnings) == 1 and "overlap" in warnings[0]


def test_phase_centres_follow_the_template_offset():
    template = single_pol("PatchV", "V", offset=(0.0, 0.1, 0.4))
    plans, _ = plan_elements(
        [element("Rx_1", "Rx", "v", 180.0, (3.0, 2.0, 0.0))],
        {"PatchV": template},
        "PatchV",
        {},
        z_ref=0.4,
    )
    assert plans[0].ports[0].phase_centre == pytest.approx((3.0, 2.0, 0.0))


def test_map_duplicates_is_exact():
    tagged = {"PortSheet__P": "PortSheet", "PortSheet1__P": "PortSheet1", "L1_Patch__P": "L1_Patch"}
    pasted = ["PortSheet1__P_1", "PortSheet__P_1", "L1_Patch__P_1"]
    assert map_duplicates(pasted, tagged) == {
        "PortSheet1__P_1": "PortSheet1",
        "PortSheet__P_1": "PortSheet",
        "L1_Patch__P_1": "L1_Patch",
    }
    with pytest.raises(RuntimeError):
        map_duplicates(["Unrelated_1"], tagged)
    with pytest.raises(RuntimeError, match="missing"):
        map_duplicates(["PortSheet__P_1"], tagged)


def test_offset_expression():
    assert offset_expression("clusterRadius", 0.0) == "clusterRadius"
    assert offset_expression("-r", -0.005) == "(-r) + (-0.0050mm)"


def test_parse_port_uses_the_last_role():
    assert parse_port("TxRxDemo1_Port_Rx_1_H") == ("Rx", "H", "TxRxDemo1_Port_Rx_1_")
    assert parse_port("HornCluster1_TxV") == ("Tx", "V", "HornCluster1_Tx")


def test_s_parameter_ports():
    assert s_parameter_ports("S(A1_Port_Rx_H,A1_Port_Tx_V)") == ("A1_Port_Rx_H", "A1_Port_Tx_V")
    assert s_parameter_ports("dB(S(a,b))") is None


def test_polarimetric_quads():
    # Hand-built horn cluster and the builder's cluster: separate H/V elements per role
    horn = ["H1_RxH", "H1_RxV", "H1_TxH", "H1_TxV"]
    assert polarimetric_quads(horn) == {
        "Rx/Tx": {"RxH": "H1_RxH", "RxV": "H1_RxV", "TxH": "H1_TxH", "TxV": "H1_TxV"}
    }
    cluster = ["C1_Port_Rx_H", "C1_Port_Rx_V", "C1_Port_Tx_H", "C1_Port_Tx_V"]
    assert list(polarimetric_quads(cluster)) == ["Rx/Tx"]
    # Dual-pol Rx elements and single-pol Tx elements with distinct labels
    ula = [f"U1_Port_Rx_{i}_{p}" for i in (1, 2) for p in "HV"] + [
        "U1_Port_Tx_1_V",
        "U1_Port_Tx_2_H",
    ]
    assert list(polarimetric_quads(ula)) == ["Rx_1/Tx", "Rx_2/Tx"]
    # Two Tx H ports and two Tx V ports without pairing: no automatic Tx group
    ula += ["U1_Port_Tx_3_V", "U1_Port_Tx_4_H"]
    assert polarimetric_quads(ula) == {}
    # Loads never count
    dual = ["D1_Port_Rx_1_H", "D1_Port_Rx_1_V", "D1_Port_Tx_1_V", "D1_Port_Tx_1_H_Load"]
    assert polarimetric_quads(dual) == {}


def test_select_quad():
    cluster = ["C1_Port_Rx_H", "C1_Port_Rx_V", "C1_Port_Tx_H", "C1_Port_Tx_V"]
    name, quad = select_quad(cluster)
    assert name == "Rx/Tx" and quad["TxV"] == "C1_Port_Tx_V"
    explicit = {"RxH": "a", "RxV": "b", "TxH": "c", "TxV": "d"}
    assert select_quad(cluster, explicit) == ("custom", explicit)
    with pytest.raises(ValueError, match="Unknown quad"):
        select_quad(cluster, "Rx_9/Tx")
    with pytest.raises(ValueError, match="no dual-polarized"):
        select_quad(["C1_Port_Rx_1_V", "C1_Port_Tx_1_V"])


def test_attach_port_table():
    templates = {"PatchV": single_pol("PatchV", "V")}
    plans, _ = plan_elements(CLUSTER, templates, "PatchV", {})
    table = sidecar_dict("PatchCluster", plans, 0.0)
    excitations = [f"PatchCluster1_{name}" for name in ("Port_Tx_H", "Port_Rx_H", "Port_Tx_V")]
    rows = attach_port_table(table, excitations)
    assert [row["excitation"] for row in rows] == excitations
    assert [row["sign"] for row in rows] == [-1, 1, -1]
    # A port the table does not know, or another design's instance, rejects the table
    assert attach_port_table(table, excitations + ["PatchCluster1_Port_Tx_9_V"]) is None
    assert attach_port_table(table, ["OtherDesign1_Port_Tx_H"]) is None
