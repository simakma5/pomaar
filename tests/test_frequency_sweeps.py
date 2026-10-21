"""Tests of the AEDT-free setup and sweep helpers shared by the builder and the SBR+ simulator."""

import pytest

from pomaar.simulator.frequency_sweeps import (
    adaptive_coverage_warning,
    adaptive_frequencies_ghz,
    design_setups,
    linked_solution_problem,
    native_sweep_arguments,
    parse_frequency_ghz,
    setup_sweeps,
    sweep_bounds_ghz,
    sweep_point_count,
    sweep_ranges,
)

# The unit cell's interpolating sweep, as PyAEDT reports its properties (with the unused
# RangeStep and RangeSamples defaults that AEDT keeps)
UNIT_CELL_SWEEP = {
    "Type": "Interpolating",
    "RangeType": "LinearCount",
    "RangeStart": "74GHz",
    "RangeEnd": "84GHz",
    "RangeCount": 81,
    "RangeStep": "1MHz",
    "RangeSamples": 11,
    "SaveFields": True,
    "SaveRadFields": False,
}


def test_parse_frequency_ghz():
    assert parse_frequency_ghz("77GHz") == 77.0
    assert parse_frequency_ghz("500MHz") == 0.5
    assert parse_frequency_ghz("79") == 79.0
    assert parse_frequency_ghz("") is None


def test_sweep_ranges_main_range_and_subranges():
    assert sweep_ranges(UNIT_CELL_SWEEP) == [
        {"RangeType": "LinearCount", "RangeStart": "74GHz", "RangeEnd": "84GHz", "RangeCount": 81}
    ]
    extra = {
        "RangeType": "LinearStep",
        "RangeStart": "84GHz",
        "RangeEnd": "86GHz",
        "RangeStep": "50MHz",
    }
    for subranges in (extra, [{**extra, "RangeCount": 401}]):
        ranges = sweep_ranges({**UNIT_CELL_SWEEP, "SweepRanges": {"Subrange": subranges}})
        assert ranges[1] == extra
        assert sweep_bounds_ghz(ranges) == (74.0, 86.0)
        assert sweep_point_count(ranges) == 81 + 41


def test_native_sweep_arguments_copy_the_ranges():
    ranges = sweep_ranges(UNIT_CELL_SWEEP)
    arguments = native_sweep_arguments("Sweep", ranges, [("Type", "Discrete")])
    assert arguments == [
        "NAME:Sweep",
        "RangeType:=", "LinearCount",
        "RangeStart:=", "74GHz",
        "RangeEnd:=", "84GHz",
        "RangeCount:=", 81,
        "Type:=", "Discrete",
    ]  # fmt: skip
    two = native_sweep_arguments("Sweep", ranges * 2)
    assert two[-1] == ["NAME:SweepRanges", ["NAME:Subrange"] + two[1:9]]
    with pytest.raises(ValueError):
        native_sweep_arguments("Sweep", [])


def test_adaptive_frequencies_of_every_solve_type():
    assert adaptive_frequencies_ghz({"SolveType": "Single", "Frequency": "79GHz"}) == [79.0]
    multi = {"SolveType": "MultiFrequency", "Frequency": "5GHz"}
    as_list = {"AdaptAt": [{"Frequency": f} for f in ("81GHz", "77GHz", "79GHz")]}
    as_keys = {"81GHz": [0.02], "77GHz": [0.02], "79GHz": [0.02]}
    for adaptive in (as_list, as_keys):
        frequencies = adaptive_frequencies_ghz({**multi, "MultipleAdaptiveFreqsSetup": adaptive})
        assert frequencies == [77.0, 79.0, 81.0]
    broadband = {
        "SolveType": "BroadBand",
        "MultipleAdaptiveFreqsSetup": {"Low": "74GHz", "High": "84GHz"},
    }
    assert adaptive_frequencies_ghz(broadband) == [74.0, 84.0]


def test_adaptive_coverage_warning_only_below_the_top_of_the_sweep():
    ranges = sweep_ranges(UNIT_CELL_SWEEP)
    inner = {"SolveType": "MultiFrequency", "MultipleAdaptiveFreqsSetup": {"77GHz": 0, "81GHz": 0}}
    assert "below the top of the sweep (84 GHz)" in adaptive_coverage_warning(inner, ranges)
    edges = {"SolveType": "MultiFrequency", "MultipleAdaptiveFreqsSetup": {"74GHz": 0, "84GHz": 0}}
    assert adaptive_coverage_warning(edges, ranges) is None
    assert adaptive_coverage_warning({}, ranges) is None


def test_only_discrete_sweeps_with_fields_can_be_linked():
    assert "adaptive" in linked_solution_problem("Setup : LastAdaptive", None)
    assert "basis frequencies" in linked_solution_problem("Setup : Sweep", UNIT_CELL_SWEEP)
    assert "'Fast'" in linked_solution_problem("Setup : Sweep", {"Type": "Fast"})
    discrete = {**UNIT_CELL_SWEEP, "Type": "Discrete", "SaveFields": False}
    assert "saves no fields" in linked_solution_problem("Setup : Sweep", discrete)
    for flags in ({"SaveFields": True}, {"SaveRadFields": "true"}):
        assert linked_solution_problem("Setup : Sweep", {**discrete, **flags}) is None


def test_saved_setups_and_sweeps_of_a_parsed_project_file():
    discrete = {**UNIT_CELL_SWEEP, "Type": "Discrete"}
    setup = {"SolveType": "MultiFrequency", "Sweeps": {"NextUniqueID": 1, "Sweep": discrete}}
    project = {
        "AnsoftProject": {
            "Created": "today",
            "HFSSModel": [
                {"Name": "Other", "AnalysisSetup": {"SolveSetups": {}}},
                {
                    "Name": "Array",
                    "AnalysisSetup": {"SolveSetups": {"NextUniqueID": 2, "Setup": setup}},
                },
            ],
        }
    }
    setups = design_setups(project, "Array")
    assert list(setups) == ["Setup"]
    # The type as saved, unlike PyAEDT's sweep objects, which always report 'Interpolating'
    assert setup_sweeps(setups["Setup"]) == {"Sweep": discrete}
    assert setup_sweeps({}) == {}
    with pytest.raises(KeyError):
        design_setups(project, "Missing")
