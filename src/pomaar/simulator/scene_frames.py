#!/usr/bin/env python3
"""
Frames of the SBR+ scene and of the processing chain (pure, no AEDT), shared by `sbr_simulator`
and `polarimetry_processor`:

* The array frame is the HFSS array design's global frame: boresight +z, H along +x, V along +y.
  The builder's port table (`<Design>_ports.json`) gives the phase centres in it.
* The scene frame is set by the boresight axis. 'x' (the default) is ISO 8855: x forward
  (boresight), y left (H), z up (V); `AntennaCS` turns the array frame onto it with the Euler
  ZYZ angles (0, 90, 90) deg, i.e. `ARRAY_TO_ISO`. 'z' keeps the array frame.
* A direction (azimuth, elevation) has its azimuth from boresight towards H and its elevation
  towards V in either frame, so the same angles steer to the same direction relative to the
  array, and frame changes leave every beamformed product unchanged.

`target_cs` holds the TargetCS definition of the scene: its spherical origin and Euler angles as
AEDT expressions in the variables targetRange, targetAzimuth, targetElevation and targetRoll.
"""

import numpy as np

BORESIGHT_AXES = ("x", "z")
DEFAULT_BORESIGHT = "x"

# Columns: the array frame's x (H), y (V) and z (boresight) axes in the ISO frame
ARRAY_TO_ISO = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
# AntennaCS Euler ZYZ angles (phi, theta, psi) in deg: Rz(0) Ry(90) Rz(90) is ARRAY_TO_ISO
ANTENNA_CS_ANGLES = {"x": ("0deg", "90deg", "90deg"), "z": ("0deg", "0deg", "0deg")}


def check_boresight(boresight):
    if boresight not in BORESIGHT_AXES:
        raise ValueError(f"boresight must be one of {BORESIGHT_AXES}, not '{boresight}'")
    return boresight


def array_to_scene(vectors, boresight=DEFAULT_BORESIGHT):
    """Array-frame vectors (..., 3) in the scene frame of `boresight` (unchanged for 'z')."""
    vectors = np.asarray(vectors, dtype=float)
    if check_boresight(boresight) == "z":
        return vectors
    return vectors @ ARRAY_TO_ISO.T


def aperture_axes(boresight=DEFAULT_BORESIGHT):
    """Indices of the scene-frame axes along H and V (the aperture plane)."""
    return (1, 2) if check_boresight(boresight) == "x" else (0, 1)


def direction_vectors(azimuth_deg, elevation_deg, boresight=DEFAULT_BORESIGHT):
    """
    Unit vectors (n_az, n_el, 3) of the direction grid in the scene frame: azimuth from boresight
    towards H, elevation towards V. With 'x', r = (cos el cos az, cos el sin az, sin el); with
    'z', r = (cos el sin az, sin el, cos el cos az).
    """
    az = np.radians(np.atleast_1d(np.asarray(azimuth_deg, dtype=float)))[:, None]
    el = np.radians(np.atleast_1d(np.asarray(elevation_deg, dtype=float)))[None, :]
    in_array_frame = np.stack(
        np.broadcast_arrays(np.cos(el) * np.sin(az), np.sin(el), np.cos(el) * np.cos(az)), axis=-1
    )
    return array_to_scene(in_array_frame, boresight)


def target_cs(boresight=DEFAULT_BORESIGHT):
    """
    TargetCS of the scene: (origin, Euler ZYZ angles, in_antenna_cs), as AEDT expressions; the
    CS is created in AntennaCS when `in_antenna_cs`, else in the global CS. Its local z points
    radially away from the antenna and, at targetRoll = 0, its local x lies along H.

    * 'x': origin targetRange (cos el cos az, cos el sin az, sin el), the processor's direction
      (az, el) (`direction_vectors`); angles (az, 90 deg - el, 90 deg + roll).
    * 'z': origin targetRange (cos az sin el, sin az sin el, cos el) in AntennaCS, with az the
      azimuth about boresight and el the polar angle off it (not the processor's angles);
      angles (az, el, roll - az), so no twist at roll = 0.
    """
    if check_boresight(boresight) == "x":
        origin = [
            "targetRange*cos(targetElevation)*cos(targetAzimuth)",
            "targetRange*cos(targetElevation)*sin(targetAzimuth)",
            "targetRange*sin(targetElevation)",
        ]
        angles = ["targetAzimuth", "90deg-targetElevation", "90deg+targetRoll"]
        return origin, angles, False
    origin = [
        "targetRange*cos(targetAzimuth)*sin(targetElevation)",
        "targetRange*sin(targetAzimuth)*sin(targetElevation)",
        "targetRange*cos(targetElevation)",
    ]
    angles = ["targetAzimuth", "targetElevation", "targetRoll - targetAzimuth"]
    return origin, angles, True


def euler_zyz(phi_deg, theta_deg, psi_deg):
    """Rotation Rz(phi) Ry(theta) Rz(psi), whose columns are the axes of an AEDT Euler ZYZ CS."""

    def rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

    def ry(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])

    phi, theta, psi = np.radians([phi_deg, theta_deg, psi_deg])
    return rz(phi) @ ry(theta) @ rz(psi)
