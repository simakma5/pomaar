"""
Analysis setups and frequency sweeps, without AEDT: shared by `hfss_array_builder` and
`sbr_simulator`, tested in tests/test_frequency_sweeps.py.

An SBR+ design receives its linked antenna's far fields only at the frequencies where the HFSS
design saved fields, and interpolates them in between. Linking an adaptive solution (fields at
the adaptive frequencies only) or an interpolating sweep (fields at its basis frequencies only)
therefore puts ripple, with slope breaks at the solved frequencies, into every SBR+ result. The
builder gives each array a discrete sweep that saves fields at every frequency, and
`sbr_simulator` links only such sweeps (`linked_solution_problem`) and copies their frequencies.
"""

import re

LINKABLE_SWEEP_TYPE = "Discrete"
RANGE_KEYS = ("RangeType", "RangeStart", "RangeEnd", "RangeCount", "RangeStep", "RangeSamples")
# AEDT keeps every range key in a sweep's properties (e.g. a default RangeStep in a LinearCount
# range); only the one that the range type uses is copied
SPACING_KEYS = {"LinearCount": "RangeCount", "LinearStep": "RangeStep", "LogScale": "RangeSamples"}
LINK_REQUIREMENT = (
    "SBR+ receives the linked antenna's far fields only at the frequencies where the HFSS "
    "design saved fields and interpolates them in between, which puts ripple with slope breaks "
    "at the solved frequencies into every result. Link a discrete sweep that saves fields "
    "(hfss_array_builder gives every array one), or add one to a hand-built design."
)


def parse_frequency_ghz(value):
    """Parses an AEDT frequency string (e.g. '77GHz', '500MHz') to GHz. Bare numbers are GHz."""
    if value is None or value == "":
        return None
    match = re.match(r"\s*([-+\d.eE]+)\s*([a-zA-Z]*)", str(value))
    if not match:
        return None
    number, unit = float(match.group(1)), match.group(2).lower()
    scale = {"": 1.0, "ghz": 1.0, "mhz": 1e-3, "khz": 1e-6, "hz": 1e-9, "thz": 1e3}
    return number * scale.get(unit, 1.0)


def _is_set(flag):
    """AEDT flags arrive as booleans or as 'true'/'false' strings."""
    return flag is True or str(flag).strip().lower() == "true"


def design_setups(project, design_name):
    """
    {setup name: properties} of the design `design_name` in a parsed project file (the dict of
    PyAEDT's `load_entire_aedt_file`), as saved: each setup's 'Sweeps' holds its sweeps.
    """
    for designs in project.get("AnsoftProject", {}).values():
        for design in designs if isinstance(designs, list) else [designs]:
            if isinstance(design, dict) and design.get("Name") == design_name:
                setups = design.get("AnalysisSetup", {}).get("SolveSetups", {})
                return {name: props for name, props in setups.items() if isinstance(props, dict)}
    raise KeyError(f"Design '{design_name}' is not in the project file.")


def setup_sweeps(setup_props):
    """{sweep name: properties} of a saved setup (see `design_setups`)."""
    sweeps = setup_props.get("Sweeps") or {}
    return {name: props for name, props in sweeps.items() if isinstance(props, dict)}


def _range(props):
    """RangeType, RangeStart, RangeEnd and the spacing key of that range type, as present."""
    keys = ("RangeType", "RangeStart", "RangeEnd", SPACING_KEYS.get(props.get("RangeType")))
    return {key: props[key] for key in keys if key and key in props}


def sweep_ranges(props):
    """
    The frequency ranges of a sweep's properties: the main range, then every 'SweepRanges'
    subrange, each a dict of RangeType, RangeStart, RangeEnd and the range's RangeCount,
    RangeStep or RangeSamples, values as AEDT gives them (e.g. '74GHz').
    """
    ranges = []
    if props.get("RangeType") and props.get("RangeStart") not in (None, ""):
        ranges.append(_range(props))
    subranges = (props.get("SweepRanges") or {}).get("Subrange", [])
    if isinstance(subranges, dict):
        subranges = [subranges]
    for subrange in subranges or []:
        ranges.append(_range(subrange))
    return ranges


def sweep_bounds_ghz(ranges):
    """(lowest, highest) frequency of `ranges` in GHz."""
    edges = [
        parse_frequency_ghz(r.get(key))
        for r in ranges
        for key in ("RangeStart", "RangeEnd")
        if r.get(key) not in (None, "")
    ]
    edges = [f for f in edges if f is not None]
    if not edges:
        raise ValueError("The sweep has no readable frequency range.")
    return min(edges), max(edges)


def sweep_point_count(ranges):
    """Number of frequencies of `ranges` (points shared by two ranges are counted twice)."""
    count = 0
    for r in ranges:
        kind = r.get("RangeType")
        if kind == "LinearCount":
            count += int(r["RangeCount"])
        elif kind == "LinearStep":
            start, stop = parse_frequency_ghz(r["RangeStart"]), parse_frequency_ghz(r["RangeEnd"])
            count += int(round((stop - start) / parse_frequency_ghz(r["RangeStep"]))) + 1
        elif kind in ("SinglePoints", "SinglePoint"):
            count += 1
        else:
            raise ValueError(f"Cannot count the frequencies of a '{kind}' range.")
    return count


def native_sweep_arguments(name, ranges, settings=()):
    """
    The argument list of AEDT's InsertFrequencySweep (and of an SBR+ setup's 'Sweeps') for a
    sweep `name` over `ranges` (see `sweep_ranges`): the first range at the top level, the
    others as 'SweepRanges' subranges, values copied as given; `settings` are further
    'Key:=' value pairs such as the sweep type.
    """
    if not ranges:
        raise ValueError(f"Sweep '{name}' needs at least one frequency range.")
    arguments = [f"NAME:{name}"]
    for key in RANGE_KEYS:
        if key in ranges[0]:
            arguments += [f"{key}:=", ranges[0][key]]
    for key, value in settings:
        arguments += [f"{key}:=", value]
    if len(ranges) > 1:
        subranges = ["NAME:SweepRanges"]
        for r in ranges[1:]:
            subrange = ["NAME:Subrange"]
            for key in RANGE_KEYS:
                if key in r:
                    subrange += [f"{key}:=", r[key]]
            subranges.append(subrange)
        arguments.append(subranges)
    return arguments


def adaptive_frequencies_ghz(setup_props):
    """
    The frequencies (GHz) at which an HFSS driven setup adapts its mesh: the single frequency,
    every multi-frequency entry, or the broadband Low and High.
    """
    solve_type = str(setup_props.get("SolveType", "")).lower()
    single = parse_frequency_ghz(setup_props.get("Frequency"))
    if solve_type == "single" and single:
        return [single]
    adaptive = setup_props.get("MultipleAdaptiveFreqsSetup", {})
    if isinstance(adaptive, dict):
        if "Low" in adaptive and "High" in adaptive:
            frequencies = [parse_frequency_ghz(adaptive[key]) for key in ("Low", "High")]
        elif isinstance(adaptive.get("AdaptAt"), list):
            frequencies = [parse_frequency_ghz(f.get("Frequency")) for f in adaptive["AdaptAt"]]
        else:
            frequencies = [parse_frequency_ghz(key) for key in adaptive]
        frequencies = sorted(f for f in frequencies if f)
        if frequencies:
            return frequencies
    return [single] if single else []


def adaptive_coverage_warning(setup_props, ranges, tolerance_ghz=1e-6):
    """
    A note when the mesh is adapted only below the top of the sweep, whose shortest wavelength
    is then not adapted for; None otherwise (or when the adaptive frequencies are unknown).
    """
    adaptive = adaptive_frequencies_ghz(setup_props)
    if not adaptive:
        return None
    top = sweep_bounds_ghz(ranges)[1]
    if max(adaptive) >= top - tolerance_ghz:
        return None
    listed = ", ".join(f"{f:g}" for f in adaptive)
    return (
        f"The mesh is adapted at {listed} GHz, below the top of the sweep ({top:g} GHz), whose "
        "shortest wavelength is then not adapted for; consider adapting at the top as well."
    )


def linked_solution_problem(solution, sweep_props):
    """
    None when `solution` ('<setup> : <sweep>', with the sweep's properties `sweep_props`, or
    None for an adaptive solution such as 'Setup : LastAdaptive') can be linked into SBR+,
    i.e. is a discrete sweep that saves fields; else a short reason why not (see
    LINK_REQUIREMENT for the full one).
    """
    if sweep_props is None:
        return "an adaptive solution, with fields only at the adaptive frequencies"
    sweep_type = sweep_props.get("Type")
    if sweep_type == "Interpolating":
        return "an interpolating sweep, with fields only at the basis frequencies it solved"
    if sweep_type != LINKABLE_SWEEP_TYPE:
        return (
            f"a '{sweep_type}' sweep; only discrete sweeps are known to give fields at every "
            "frequency"
        )
    if not (_is_set(sweep_props.get("SaveFields")) or _is_set(sweep_props.get("SaveRadFields"))):
        return "a discrete sweep that saves no fields"
    return None
