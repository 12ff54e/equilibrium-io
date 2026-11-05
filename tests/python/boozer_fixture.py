"""Manufactured Boozer container checks, independent of plotting."""

import struct

import numpy as np

from equilibrium_io.boozer import (
    ASYMMETRIC_FAMILY_NAMES, ASYMMETRIC_FOURIER_CONVENTION,
    COORDINATE_CONVENTION, FAMILY_NAMES, FOURIER_CONVENTION, PERIOD,
    MAGNETIC_CONVENTION, PROFILE_CONVENTION, TRACING_CONVENTION,
)


def _string(value):
    encoded = value.encode("utf-8")
    return struct.pack("<i", len(encoded)) + encoded


def _manufactured_file(path, version=3, lasym=False):
    source_ns = 4
    first_surface = 1
    surfaces = source_ns - first_surface
    ntheta, nzeta = 12, 6
    mode_m = np.array([0, 1, 0, 1], dtype="<i4")
    mode_n = np.array([0, 0, 1, 1], dtype="<i4")
    s = np.array([1.0 / 3.0, 2.0 / 3.0, 1.0], dtype="<f8")
    iota = np.array([0.0, 0.45, 0.60], dtype="<f8")
    radius = np.array([0.15, 0.30, 0.45])
    names = FAMILY_NAMES + (ASYMMETRIC_FAMILY_NAMES if lasym else ())
    families = {name: np.zeros((surfaces, 4), dtype="<f8") for name in names}
    families["rmncc"][:, 0] = 1.5
    families["rmncc"][:, 1] = radius
    families["zmnsc"][:, 1] = radius
    families["numnsc"][:, 1] = 0.08
    families["numncs"][:, 2] = 0.03
    if lasym:
        families["rmnsc"][:, 1] = 0.04
        families["zmncc"][:, 0] = 0.11
    theta = PERIOD * np.arange(ntheta) / ntheta
    alpha = PERIOD * np.arange(nzeta) / nzeta
    b = np.empty((surfaces, nzeta, ntheta), dtype="<f8")
    for surface in range(surfaces):
        b[surface] = (
            1.0 + 0.1 * surface + 0.2 * np.cos(theta)[None, :]
            + 0.05 * np.cos(alpha)[:, None]
        )
    sqrtg = -1.0 / (b * b)
    b2j00 = -np.ones(surfaces, dtype="<f8")

    header = [
        10 if lasym else 8, source_ns, ntheta, nzeta, 2, 0, 5, first_surface,
        ntheta, nzeta, 1, 1, 2 if version >= 5 else 4,
    ]
    payload = bytearray(f"MCBOOZ0{version}".encode())
    payload += struct.pack("<i", version)
    payload += _string(COORDINATE_CONVENTION)
    payload += _string(ASYMMETRIC_FOURIER_CONVENTION if lasym
                       else FOURIER_CONVENTION)
    payload += _string("manufactured.json")
    payload += struct.pack("<13i", *header)
    if version >= 5:
        payload += struct.pack("<i", int(lasym))
    payload += struct.pack("<d", 1.0e-12)
    for array in (s, iota, mode_m, mode_n):
        payload += array.tobytes()
    for name in names:
        payload += families[name].tobytes()
    for array in (b, sqrtg, b2j00):
        payload += array.tobytes()
    if version >= 5:
        payload += _string(TRACING_CONVENTION) + struct.pack("<id", 1, 2.0 * np.pi)
        phip = -np.ones(source_ns)
        chip = phip * np.r_[0.3, iota]
        for array in ([1.0], phip, chip, [1.5], [0.0],
                      [0.11 if lasym else 0.0], [0.0], [1/6, 0.5]):
            payload += np.asarray(array, dtype="<f8").tobytes()
        for field in range(7):
            payload += np.full((2, nzeta, ntheta), field + 1.0,
                               dtype="<f8").tobytes()
    if version >= 6:
        payload += _string(MAGNETIC_CONVENTION)
        for field in range(7):
            payload += np.full((surfaces, nzeta, ntheta), field+1.0,
                               dtype="<f8").tobytes()
        for displacement in (0.13, 0.07):
            payload += np.full_like(b, displacement).tobytes()
    if version >= 7:
        payload += _string(PROFILE_CONVENTION)
        for values in ([125000.0, 62500.0, 0.0], [1.0, 2.0, 3.0], [-5.0, -6.0, -7.0]):
            payload += np.asarray(values, dtype="<f8").tobytes()
    if version >= 8:
        payload += _string("list")
    path.write_bytes(payload)
    return radius, b


