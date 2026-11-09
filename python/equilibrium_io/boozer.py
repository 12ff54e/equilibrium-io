"""Read the geometry and field-strength subset of Boozer v7/v8 files.

This is an analysis reader, not a lossless record reader. It returns Fourier
geometry, iota and |B|. It does not retain optional tracing, magnetic component
or surface-profile records. Use the C++ API for complete records.
"""

from dataclasses import dataclass
import struct

import numpy as np


MAGICS = {f"MCBOOZ0{version}".encode(): version for version in range(7, 9)}
PERIOD = 2.0 * np.pi
FAMILY_NAMES = (
    "rmncc", "rmnss", "zmnsc", "zmncs", "numnsc", "numncs",
)
ASYMMETRIC_FAMILY_NAMES = (
    "rmnsc", "rmncs", "zmncc", "zmnss", "numncc", "numnss",
)
SCHEMAS = {f"magnetic-coordinate-boozer-v{version}": version
           for version in range(7, 9)}
COORDINATE_CONVENTION = (
    "mixed-grid-v2: theta_b uniform; alpha=nfp*zeta is the unchanged source "
    "field-period angle; alpha_b=alpha+nfp*nu; nu is in physical toroidal "
    "radians"
)
FOURIER_CONVENTION = (
    "real-parity-v3: f=sum[cc*cos(m*theta_b)*cos(n*alpha) + "
    "ss*sin(m*theta_b)*sin(n*alpha)] for even fields and "
    "f=sum[sc*sin(m*theta_b)*cos(n*alpha) + "
    "cs*cos(m*theta_b)*sin(n*alpha)] for odd fields; m,n are nonnegative "
    "and alpha is a field-period angle"
)
ASYMMETRIC_FOURIER_CONVENTION = (
    "real-parity-v4: f=sum[cc*cos(m*theta_b)*cos(n*alpha) + "
    "ss*sin(m*theta_b)*sin(n*alpha) + sc*sin(m*theta_b)*cos(n*alpha) + "
    "cs*cos(m*theta_b)*sin(n*alpha)]; m,n are nonnegative "
    "and alpha is a field-period angle; no stellarator symmetry is imposed"
)
TRACING_CONVENTION = (
    "tracing-v1: signed flux/(2*pi) derivatives in source s; axis m=0 "
    "geometry in alpha=nfp*zeta; axis_half fields are unmodified native "
    "half-grid samples in (s,theta,zeta), not Boozer axis values"
)

MAGNETIC_CONVENTION = (
    "boozer-magnetic-v1: boozer fields use (s,theta_b,zeta_b) with physical "
    "zeta_b; sampled on the mixed grid with theta=theta_b+source_theta_shift "
    "and zeta_b=alpha/nfp+nu_grid; shifts are unwrapped and untruncated"
)

PROFILE_CONVENTION = (
    "surface-profiles-v1: pressure in Pa; B=beta*grad(s)+I*grad(theta_b)+G*grad(zeta_b); "
    "current_i=I and current_g=G in T m with physical zeta_b; all profiles share the non-axis s/iota grid"
)
PROFILE_UNITS = {"pressure": "Pa", "current_i": "T m", "current_g": "T m"}


def _validate_profiles(profiles, iota, path):
    for values in profiles.values():
        if values.shape != iota.shape or not np.all(np.isfinite(values)):
            raise ValueError(f"error: invalid Boozer surface profiles in {path}")
    iota_i = iota * profiles["current_i"]
    g = profiles["current_g"]
    normalization = g + iota_i
    if np.any(~np.isfinite(normalization)) or np.any(
            np.abs(normalization) <= 64 * np.finfo(float).eps *
            (np.abs(g) + np.abs(iota_i))):
        raise ValueError(f"error: ill-conditioned Boozer surface profiles in {path}")


def _family_names(lasym):
    if lasym not in (0, 1):
        raise ValueError("error: invalid Boozer symmetry flag")
    return FAMILY_NAMES + (ASYMMETRIC_FAMILY_NAMES if lasym else ())


@dataclass(frozen=True)
class BoozerData:
    """Validated Boozer geometry and field-strength subset for analysis."""

    source_path: str
    source_ns: int
    nfp: int
    first_surface: int
    ntheta: int
    nzeta: int
    s: np.ndarray
    iota: np.ndarray
    mode_m: np.ndarray
    mode_n: np.ndarray
    families: dict
    b: np.ndarray

    @property
    def surface_count(self):
        return self.source_ns - self.first_surface

    @property
    def theta_b(self):
        return PERIOD * np.arange(self.ntheta) / self.ntheta

    @property
    def alpha(self):
        return PERIOD * np.arange(self.nzeta) / self.nzeta

    @property
    def zeta(self):
        """Alias for the stored field-period angle, retained for callers."""
        return self.alpha


class _BinaryReader:
    def __init__(self, path):
        try:
            self._stream = open(path, "rb")
        except OSError as exc:
            raise OSError(f"could not read Boozer file {path}: {exc}") \
                from exc
        self.path = path
        self._stream.seek(0, 2)
        self._size = self._stream.tell()
        self._stream.seek(0)

    def close(self):
        self._stream.close()

    def read(self, size):
        if size < 0 or size > self._size - self._stream.tell():
            raise ValueError(f"error: truncated Boozer file {self.path}")
        value = self._stream.read(size)
        if len(value) != size:
            raise ValueError(f"error: truncated Boozer file {self.path}")
        return value

    def i32(self):
        return struct.unpack("<i", self.read(4))[0]

    def f64(self):
        return struct.unpack("<d", self.read(8))[0]

    def string(self):
        size = self.i32()
        if not 0 <= size < 1 << 24:
            raise ValueError(
                f"error: corrupt string length in Boozer file {self.path}")
        try:
            return self.read(size).decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"error: invalid UTF-8 in Boozer file {self.path}") from exc

    def array(self, dtype, count):
        itemsize = np.dtype(dtype).itemsize
        return np.frombuffer(self.read(itemsize * count), dtype=dtype).copy()


def _load_binary(path):
    reader = _BinaryReader(path)
    try:
        version = MAGICS.get(reader.read(8))
        if version is None or reader.i32() != version:
            raise ValueError(f"error: unsupported Boozer binary format in {path}")
        coordinate_convention = reader.string()
        fourier_convention = reader.string()
        if coordinate_convention != COORDINATE_CONVENTION:
            raise ValueError(
                f"error: unsupported Boozer coordinate convention in {path}")
        source_path = reader.string()
        header = [reader.i32() for _ in range(13)]
        (source_format_version, source_ns, source_ntheta, source_nzeta,
         source_mpol, source_ntor, nfp, first_surface, ntheta, nzeta,
         mmax, nmax, radial_order) = header
        lasym = reader.i32()
        family_names = _family_names(lasym)
        expected_fourier = (ASYMMETRIC_FOURIER_CONVENTION if lasym
                            else FOURIER_CONVENTION)
        if fourier_convention != expected_fourier:
            raise ValueError(f"error: unsupported Boozer Fourier convention in {path}")
        resonance_tolerance = reader.f64()
        del source_format_version, source_mpol, resonance_tolerance
        if mmax < 0 or nmax < 0:
            raise ValueError(f"error: invalid Boozer mode bounds in {path}")
        surfaces = source_ns - first_surface
        modes = (mmax + 1) * (nmax + 1)
        if not (1 <= surfaces <= 1 << 20 and 1 <= modes <= 1 << 20 and
                ntheta >= 1 and nzeta >= 1):
            raise ValueError(f"error: invalid Boozer dimensions in {path}")
        s = reader.array("<f8", surfaces)
        iota = reader.array("<f8", surfaces)
        mode_m = reader.array("<i4", modes)
        mode_n = reader.array("<i4", modes)
        families = {
            name: reader.array("<f8", surfaces * modes).reshape(surfaces, modes)
            for name in family_names
        }
        b = reader.array("<f8", surfaces * nzeta * ntheta).reshape(
            surfaces, nzeta, ntheta)
        reader.array("<f8", surfaces * nzeta * ntheta)  # sqrt_g_b
        reader.array("<f8", surfaces)  # b2j00
        if reader.string() != TRACING_CONVENTION:
            raise ValueError(f"error: unsupported Boozer tracing convention in {path}")
        aphi_count = reader.i32()
        if not (1 <= aphi_count <= 1 << 20 and radial_order in (2, 4)
                and radial_order < source_ns <= 1 << 20
                and 0 <= source_ntor <= 1 << 15
                and source_ntheta > 0 and source_nzeta > 0):
            raise ValueError(f"error: invalid Boozer tracing dimensions in {path}")
        # Plotting uses the non-axis grid. Still check that the required
        # tracing trailer is complete; field reconstruction uses the library API.
        count = (1 + aphi_count + 2 * source_ns + 4 * (source_ntor + 1)
                 + radial_order
                 + 7 * radial_order * source_ntheta * source_nzeta)
        if not np.all(np.isfinite(reader.array("<f8", count))):
            raise ValueError(f"error: non-finite Boozer tracing data in {path}")
        if reader.string() != MAGNETIC_CONVENTION:
            raise ValueError(f"error: unsupported Boozer magnetic convention in {path}")
        count = 9 * surfaces * ntheta * nzeta
        if not np.all(np.isfinite(reader.array("<f8", count))):
            raise ValueError(f"error: non-finite Boozer magnetic data in {path}")
        if reader.string() != PROFILE_CONVENTION:
            raise ValueError(f"error: unsupported Boozer profile convention in {path}")
        profiles = {name: reader.array("<f8", surfaces)
                    for name in PROFILE_UNITS}
        _validate_profiles(profiles, iota, path)
        if version == 8 and reader.string() not in ("power_series", "list"):
            raise ValueError(f"error: unsupported Boozer flux profile in {path}")
    finally:
        reader.close()
    return BoozerData(source_path, source_ns, nfp, first_surface, ntheta,
                      nzeta, s, iota, mode_m, mode_n, families, b)


def _decode(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _load_netcdf(path):
    try:
        from scipy.io import netcdf_file
        with netcdf_file(path, "r", mmap=False) as nc:
            version = SCHEMAS.get(_decode(getattr(nc, "schema", "")))
            if version is None:
                raise ValueError(f"error: unsupported Boozer NetCDF schema in {path}")

            def attr(name):
                return int(getattr(nc, name))

            lasym = attr("lasym")
            source_path = _decode(getattr(nc, "source_path", ""))
            source_ns = attr("source_ns")
            nfp = attr("nfp")
            first_surface = attr("first_surface")
            ntheta = attr("ntheta")
            nzeta = attr("nzeta")
            arrays = {name: np.asarray(nc.variables[name][:]).copy()
                      for name in ("s", "iota", "mode_m", "mode_n", "B")}
            families = {name: np.asarray(nc.variables[name][:], dtype=float).copy()
                        for name in _family_names(lasym)}
            if _decode(getattr(nc, "profile_convention", "")) != PROFILE_CONVENTION:
                raise ValueError(f"error: unsupported Boozer profile convention in {path}")
            profiles = {}
            for name, unit in PROFILE_UNITS.items():
                variable = nc.variables[name]
                if _decode(getattr(variable, "units", "")) != unit:
                    raise ValueError(f"error: invalid Boozer profile units in {path}")
                profiles[name] = np.asarray(variable[:], dtype=float).copy()
            _validate_profiles(profiles, arrays["iota"], path)
    except OSError as exc:
        raise OSError(f"could not read Boozer file {path}: {exc}") \
            from exc
    return BoozerData(source_path, source_ns, nfp, first_surface, ntheta,
                      nzeta, np.asarray(arrays["s"], dtype=float),
                      np.asarray(arrays["iota"], dtype=float),
                      np.asarray(arrays["mode_m"], dtype=np.int32),
                      np.asarray(arrays["mode_n"], dtype=np.int32), families,
                      np.asarray(arrays["B"], dtype=float))


def _load_hdf5(path):
    try:
        import h5py
        with h5py.File(path, "r") as h5:
            version = SCHEMAS.get(_decode(h5.attrs.get("schema", "")))
            if version is None:
                raise ValueError(f"error: unsupported Boozer HDF5 schema in {path}")
            lasym = int(h5.attrs["lasym"])
            source_path = _decode(h5.attrs.get("source_path", ""))
            source_ns = int(h5.attrs["source_ns"])
            nfp = int(h5.attrs["nfp"])
            first_surface = int(h5.attrs["first_surface"])
            ntheta = int(h5.attrs["ntheta"])
            nzeta = int(h5.attrs["nzeta"])
            s = np.asarray(h5["s"][:], dtype=float)
            iota = np.asarray(h5["iota"][:], dtype=float)
            mode_m = np.asarray(h5["mode_m"][:], dtype=np.int32)
            mode_n = np.asarray(h5["mode_n"][:], dtype=np.int32)
            families = {name: np.asarray(h5[name][:], dtype=float)
                        for name in _family_names(lasym)}
            b = np.asarray(h5["B"][:], dtype=float)
            if _decode(h5.attrs.get("profile_convention", "")) != PROFILE_CONVENTION:
                raise ValueError(f"error: unsupported Boozer profile convention in {path}")
            profiles = {}
            for name, unit in PROFILE_UNITS.items():
                variable = h5[name]
                if _decode(variable.attrs.get("units", "")) != unit:
                    raise ValueError(f"error: invalid Boozer profile units in {path}")
                profiles[name] = np.asarray(variable[:], dtype=float)
            _validate_profiles(profiles, iota, path)
    except OSError as exc:
        raise OSError(f"could not read Boozer file {path}: {exc}") \
            from exc
    return BoozerData(source_path, source_ns, nfp, first_surface, ntheta,
                      nzeta, s, iota, mode_m, mode_n, families, b)


def _validate(data, path):
    surfaces = data.source_ns - data.first_surface
    modes = data.mode_m.size
    if data.nfp < 1 or data.first_surface < 1 or surfaces < 1 or \
            data.ntheta < 1 or data.nzeta < 1:
        raise ValueError(f"error: invalid Boozer dimensions in {path}")
    if data.s.shape != (surfaces,) or data.iota.shape != (surfaces,) or \
            data.mode_n.shape != (modes,) or data.b.shape != (
                surfaces, data.nzeta, data.ntheta):
        raise ValueError(f"error: inconsistent Boozer array extents in {path}")
    if modes < 1 or any(np.asarray(data.families[name]).shape != (surfaces, modes)
                        for name in data.families):
        raise ValueError(f"error: inconsistent Boozer spectrum extents in {path}")
    real_arrays = [data.s, data.iota, data.b]
    real_arrays.extend(data.families.values())
    if any(not np.all(np.isfinite(array)) for array in real_arrays):
        raise ValueError(f"error: non-finite data in Boozer file {path}")
    if np.any(data.b <= 0.0) or np.any(np.diff(data.s) <= 0.0):
        raise ValueError(f"error: nonphysical Boozer grid in {path}")
    if np.any(data.mode_m < 0) or np.any(data.mode_n < 0):
        raise ValueError(f"error: negative Boozer mode number in {path}")
    return data


def load_boozer(path):
    """Read a geometry/|B| subset from a Boozer v7/v8 container.

    Binary input needs NumPy; NetCDF and HDF5 input also need scipy and h5py,
    respectively. Malformed data raises ValueError; file errors raise OSError.
    Optional record groups are not returned. This function cannot round-trip
    complete Boozer files.
    """
    try:
        with open(path, "rb") as stream:
            magic = stream.read(8)
    except OSError as exc:
        raise OSError(f"could not read Boozer file {path}: {exc}") \
            from exc
    try:
        if magic in MAGICS:
            data = _load_binary(path)
        elif magic.startswith(b"CDF"):
            data = _load_netcdf(path)
        elif magic.startswith(b"\x89HDF"):
            data = _load_hdf5(path)
        else:
            raise ValueError(f"error: unrecognized Boozer container {path}")
        return _validate(data, path)
    except (KeyError, AttributeError, TypeError, struct.error) as error:
        raise ValueError(f"invalid Boozer record in {path}: {error}") from error


