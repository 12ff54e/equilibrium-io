"""Boozer analysis readers own their format and malformed-input checks."""

import importlib.util
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from equilibrium_io import load_boozer
from equilibrium_io.boozer import PROFILE_CONVENTION, PROFILE_UNITS
from boozer_fixture import _manufactured_file


class BoozerIoTest(unittest.TestCase):
    def test_binary_versions_and_symmetry(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "boozer.bin"
            for version, lasym in ((3, False), (4, True), (5, False), (5, True),
                                   (6, False), (6, True), (7, False), (7, True),
                                   (8, False), (8, True)):
                with self.subTest(version=version, lasym=lasym):
                    radius, expected_b = _manufactured_file(path, version, lasym)
                    data = load_boozer(path)
                    self.assertEqual((data.source_ns, data.nfp), (4, 5))
                    self.assertEqual(data.source_path, "manufactured.json")
                    self.assertEqual(len(data.families), 12 if lasym else 6)
                    np.testing.assert_array_equal(data.b, expected_b)
                    np.testing.assert_array_equal(data.families["rmncc"][:, 1], radius)
                    np.testing.assert_array_equal(data.iota, [0, 0.45, 0.60])
                    path.write_bytes(path.read_bytes()[:-1])
                    with self.assertRaisesRegex(ValueError, "truncated"):
                        load_boozer(path)

    def test_profile_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "boozer.bin"
            _manufactured_file(path, 7)
            original = path.read_bytes()
            for invalid in ("convention", "nonfinite", "normalization"):
                payload = bytearray(original)
                if invalid == "convention":
                    payload[payload.index(PROFILE_CONVENTION.encode())] = ord("x")
                elif invalid == "nonfinite":
                    struct.pack_into("<d", payload, len(payload) - 9 * 8, float("nan"))
                else:
                    struct.pack_into("<d", payload, len(payload) - 3 * 8, 0.0)
                path.write_bytes(payload)
                with self.subTest(invalid=invalid), self.assertRaisesRegex(
                        ValueError, "profile"):
                    load_boozer(path)

    def test_file_errors_do_not_exit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "missing.bin"
            with self.assertRaises(OSError):
                load_boozer(path)
            path.write_bytes(b"not-boozer")
            with self.assertRaisesRegex(ValueError, "unrecognized"):
                load_boozer(path)

    def test_huge_extent_rejected_before_allocation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "huge.bin"
            _manufactured_file(path)
            payload = bytearray(path.read_bytes())
            offset = 12
            for _ in range(3):
                size, = struct.unpack_from("<i", payload, offset)
                offset += 4 + size
            struct.pack_into("<i", payload, offset + 8 * 4, 2**30)
            path.write_bytes(payload)
            with self.assertRaisesRegex(ValueError, "truncated"):
                load_boozer(path)

    def check_codec(self, codec):
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / "boozer.bin"
            _manufactured_file(binary, 7, True)
            expected = load_boozer(binary)
            path = Path(directory) / ("boozer." + codec)
            attrs = dict(schema="magnetic-coordinate-boozer-v7", lasym=1,
                         source_path="manufactured.json", source_ns=4,
                         nfp=5, first_surface=1, ntheta=12, nzeta=6,
                         profile_convention=PROFILE_CONVENTION)
            arrays = {"s": (expected.s, ("s",)), "iota": (expected.iota, ("s",)),
                      "mode_m": (expected.mode_m, ("mode",)),
                      "mode_n": (expected.mode_n, ("mode",)),
                      "B": (expected.b, ("s", "zeta", "theta"))}
            arrays.update({name: (values, ("s", "mode"))
                           for name, values in expected.families.items()})
            arrays.update({name: (np.asarray(values), ("s",)) for name, values in (
                ("pressure", [125000., 62500., 0.]), ("current_i", [1., 2., 3.]),
                ("current_g", [-5., -6., -7.]))})
            if codec == "nc":
                from scipy.io import netcdf_file
                with netcdf_file(path, "w") as output:
                    for name, count in (("s", 3), ("mode", 4), ("theta", 12), ("zeta", 6)):
                        output.createDimension(name, count)
                    for name, value in attrs.items():
                        setattr(output, name, value)
                    for name, (values, dimensions) in arrays.items():
                        variable = output.createVariable(name, "i" if name.startswith("mode_") else "d", dimensions)
                        variable[:] = values
                        if name in PROFILE_UNITS:
                            variable.units = PROFILE_UNITS[name]
            else:
                import h5py
                with h5py.File(path, "w") as output:
                    output.attrs.update(attrs)
                    for name, (values, _) in arrays.items():
                        variable = output.create_dataset(name, data=values)
                        if name in PROFILE_UNITS:
                            variable.attrs["units"] = PROFILE_UNITS[name]
            actual = load_boozer(path)
            np.testing.assert_array_equal(actual.b, expected.b)
            for name in expected.families:
                np.testing.assert_array_equal(actual.families[name], expected.families[name])
            if codec == "h5":
                with h5py.File(path, "r+") as output:
                    del output["mode_n"]
                with self.assertRaisesRegex(ValueError, "invalid Boozer record"):
                    load_boozer(path)

    @unittest.skipUnless(importlib.util.find_spec("scipy"), "NetCDF extra not installed")
    def test_netcdf(self):
        self.check_codec("nc")

    @unittest.skipUnless(importlib.util.find_spec("h5py"), "HDF5 extra not installed")
    def test_hdf5(self):
        self.check_codec("h5")


if __name__ == "__main__":
    unittest.main()
