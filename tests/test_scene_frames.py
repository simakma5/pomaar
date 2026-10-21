"""Tests of the scene frames shared by the SBR+ scene and the processor (no AEDT needed)."""

import math
import re

import numpy as np
import pytest

from pomaar.simulator.scene_frames import (
    ANTENNA_CS_ANGLES,
    ARRAY_TO_ISO,
    array_to_scene,
    direction_vectors,
    euler_zyz,
    target_cs,
    target_steering_angles,
)

H, V, BORESIGHT = np.eye(3)  # array-frame axes


def evaluate(expression, variables):
    """Value of an AEDT expression of the target variables (angles in deg, 'deg' units)."""
    python = re.sub(r"(\d+(?:\.\d+)?)deg", r"(\1*pi/180)", expression)
    names = {"cos": math.cos, "sin": math.sin, "pi": math.pi}
    return eval(python, names, variables)  # noqa: S307 (fixed expressions of this module)


def target_frame(boresight, azimuth, elevation, roll=0.0, target_range=4.0):
    """TargetCS origin and axes (columns) in the scene frame, from `target_cs`'s expressions."""
    variables = {
        "targetRange": target_range,
        "targetAzimuth": math.radians(azimuth),
        "targetElevation": math.radians(elevation),
        "targetRoll": math.radians(roll),
    }
    origin_expressions, angle_expressions, in_antenna_cs = target_cs(boresight)
    origin = np.array([evaluate(e, variables) for e in origin_expressions])
    axes = euler_zyz(*np.degrees([evaluate(e, variables) for e in angle_expressions]))
    if in_antenna_cs:
        antenna = euler_zyz(*[float(a.removesuffix("deg")) for a in ANTENNA_CS_ANGLES[boresight]])
        origin, axes = antenna @ origin, antenna @ axes
    return origin, axes


def test_antenna_cs_turns_the_array_frame_onto_the_scene():
    angles = [float(a.removesuffix("deg")) for a in ANTENNA_CS_ANGLES["x"]]
    np.testing.assert_allclose(euler_zyz(*angles), ARRAY_TO_ISO, atol=1e-12)
    # H -> +y (left), V -> +z (up), boresight -> +x (forward)
    np.testing.assert_allclose(array_to_scene([H, V, BORESIGHT]), np.eye(3)[[1, 2, 0]])
    np.testing.assert_allclose(array_to_scene([H, V, BORESIGHT], "z"), np.eye(3))


@pytest.mark.parametrize("boresight", ["x", "z"])
@pytest.mark.parametrize("azimuth, elevation", [(0, 0), (20, 0), (-35, 10), (50, -25)])
def test_target_cs_faces_the_antenna(boresight, azimuth, elevation):
    origin, axes = target_frame(boresight, azimuth, elevation)
    assert np.linalg.norm(origin) == pytest.approx(4.0)
    np.testing.assert_allclose(axes.T @ axes, np.eye(3), atol=1e-12)
    # Local z radially outward: the targets face the antenna along local -z
    np.testing.assert_allclose(axes[:, 2], origin / 4.0, atol=1e-12)


@pytest.mark.parametrize("azimuth, elevation", [(0, 0), (20, 0), (-35, 10), (50, -25)])
def test_x_boresight_target_sits_where_the_processor_steers(azimuth, elevation):
    origin, axes = target_frame("x", azimuth, elevation)
    np.testing.assert_allclose(origin / 4.0, direction_vectors(azimuth, elevation)[0, 0])
    # At zero roll, local x is horizontal (along H at boresight, turned with the azimuth)
    assert axes[2, 0] == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_allclose(
        axes[:, 0],
        [-math.sin(math.radians(azimuth)), math.cos(math.radians(azimuth)), 0.0],
        atol=1e-12,
    )


@pytest.mark.parametrize("boresight", ["x", "z"])
def test_target_cs_at_boresight_is_aligned_with_h_and_v(boresight):
    origin, axes = target_frame(boresight, 0.0, 0.0)
    scene = array_to_scene([H, V, BORESIGHT], boresight)
    np.testing.assert_allclose(origin, 4.0 * scene[2], atol=1e-12)
    np.testing.assert_allclose(axes, scene.T, atol=1e-12)


@pytest.mark.parametrize("boresight", ["x", "z"])
def test_target_roll_turns_about_the_line_of_sight(boresight):
    origin, axes = target_frame(boresight, 15.0, 8.0)
    rolled_origin, rolled = target_frame(boresight, 15.0, 8.0, roll=30.0)
    np.testing.assert_allclose(rolled_origin, origin)
    np.testing.assert_allclose(rolled[:, 2], axes[:, 2], atol=1e-12)
    cosine = rolled[:, 0] @ axes[:, 0]
    assert cosine == pytest.approx(math.cos(math.radians(30.0)))


def test_directions_agree_between_frames():
    azimuth, elevation = np.arange(-60, 61, 15.0), np.arange(-30, 31, 10.0)
    np.testing.assert_allclose(
        direction_vectors(azimuth, elevation, "x"),
        array_to_scene(direction_vectors(azimuth, elevation, "z")),
    )
    # Azimuth towards H and elevation towards V in both frames
    for boresight in ("x", "z"):
        scene = array_to_scene([H, V, BORESIGHT], boresight)
        np.testing.assert_allclose(
            direction_vectors(90.0, 0.0, boresight)[0, 0], scene[0], atol=1e-12
        )
        np.testing.assert_allclose(
            direction_vectors(0.0, 90.0, boresight)[0, 0], scene[1], atol=1e-12
        )


@pytest.mark.parametrize("boresight", ["x", "z"])
@pytest.mark.parametrize("azimuth, elevation", [(0, 0), (20, 5), (-35, 10), (90, 10), (150, 30)])
def test_steering_angles_point_at_the_target(boresight, azimuth, elevation):
    # The processor steered to these angles looks exactly at TargetCS's origin
    origin, _ = target_frame(boresight, azimuth, elevation)
    steer = target_steering_angles(azimuth, elevation, boresight)
    np.testing.assert_allclose(direction_vectors(*steer, boresight)[0, 0], origin / 4.0, atol=1e-12)


def test_z_scene_angles_are_spherical():
    assert target_steering_angles(90.0, 10.0, "z") == pytest.approx((0.0, 10.0))
    assert target_steering_angles(0.0, 10.0, "z") == pytest.approx((10.0, 0.0))
    assert target_steering_angles(26.7402, 11.169, "z") == pytest.approx((10.0, 5.0), abs=1e-3)
