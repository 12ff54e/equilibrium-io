#!/usr/bin/env python3
"""Read C++ typed-input fixtures through the independent Python package.

Usage: python3 tests/python/test_state_io.py /path/to/equilibrium-io/build/tests/equilibrium_asymmetric_io
"""

from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from equilibrium_io import ASYM_FAM_NAMES, FAM_NAMES, load_state


class StateIoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="equilibrium_io_")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        result = subprocess.run([str(FIXTURE_WRITER), str(cls.directory)],
                                capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)

    def test_all_containers_preserve_large_and_empty_inputs(self):
        for name in ("large", "empty"):
            for suffix in (".bin", ".ckpt", ".nc", ".h5"):
                path = self.directory / (name + suffix)
                if not path.exists():  # optional C++ output backend
                    continue
                with self.subTest(name=name, suffix=suffix):
                    ns, mnmax, families, params, _ = load_state(path)
                    self.assertEqual((ns, mnmax), (5, 6))
                    self.assertEqual(tuple(families), FAM_NAMES + ASYM_FAM_NAMES)
                    for index, values in enumerate(families.values()):
                        np.testing.assert_array_equal(
                            values, index * 100.0 + np.arange(30) + 0.125)
                    self.assertTrue(params["lasym"])
                    self.assertEqual(params["raxis_s"], [0.0, 0.002])
                    self.assertEqual(params["zaxis_c"], [0.03, -0.001])
                    rbs = [] if name == "empty" else [
                        (1, 1, 0.01), (2, -1, 0.002), (1, 1, -0.003),
                    ] + [(2, -1, 0.0)] * 10000
                    zbc = [] if name == "empty" else [
                        (0, 0, 0.04), (1, -1, 0.03), (0, 1, 0.001),
                    ] + [(1, 1, 0.0)] * 10000
                    self.assertEqual(params["rbs"], rbs)
                    self.assertEqual(params["zbc"], zbc)

    def test_previous_versions_and_symmetry(self):
        for name, lasym in (("large", True), ("symmetric", False)):
            for suffix, version in ((".bin", 10), (".ckpt", 8)):
                source = self.directory / (name + suffix)
                payload = bytearray(source.read_bytes()[:-16])
                if not lasym:
                    payload += bytes(13 * 4)
                struct.pack_into("<i", payload, 8, version)
                path = self.directory / ("previous-" + name + suffix)
                path.write_bytes(payload)
                with self.subTest(name=name, version=version):
                    data = load_state(path)
                    self.assertEqual(len(data[2]), 12 if lasym else 6)
                    self.assertEqual(data[3]["lasym"], lasym)
                    self.assertEqual(data[3]["pflux_type"], "power_series")
                    self.assertEqual(data[3]["rbc"], load_state(source)[3]["rbc"])

    def test_flux_selector_extension(self):
        for suffix, version in ((".bin", 11), (".ckpt", 9)):
            payload = bytearray((self.directory / ("large" + suffix)).read_bytes())
            struct.pack_into("<i", payload, 8, version)
            payload = payload[:-16] + struct.pack("<i", 4) + b"list"
            path = self.directory / ("flux-list" + suffix)
            path.write_bytes(payload)
            self.assertEqual(load_state(path)[3]["pflux_type"], "list")
            path.write_bytes(payload[:-1])
            with self.assertRaisesRegex(ValueError, "truncated string"):
                load_state(path)

    def test_json_record_versions_are_rejected(self):
        for suffix, version in ((".bin", 9), (".ckpt", 7)):
            path = self.directory / ("old" + suffix)
            payload = bytearray((self.directory / ("large" + suffix)).read_bytes())
            struct.pack_into("<i", payload, 8, version)
            path.write_bytes(payload)
            with self.subTest(suffix=suffix), self.assertRaisesRegex(
                    ValueError, "unsupported"):
                load_state(path)

    def test_partial_hdf5_boundary_is_rejected(self):
        source = self.directory / "empty.h5"
        if not source.exists():
            self.skipTest("HDF5 output backend unavailable")
        import h5py
        path = self.directory / "partial.h5"
        shutil.copyfile(source, path)
        with h5py.File(path, "r+") as output:
            del output["zbc_n"]
        with self.assertRaisesRegex(ValueError, "incomplete asymmetric boundary"):
            load_state(path)

    def test_truncated_input_is_rejected(self):
        source = self.directory / "large.ckpt"
        path = self.directory / "truncated.ckpt"
        path.write_bytes(source.read_bytes()[:-8])
        with self.assertRaises((ValueError, struct.error)):
            load_state(path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError(__doc__)
    FIXTURE_WRITER = Path(sys.argv.pop(1)).resolve()
    unittest.main()
