# Boozer output v8: flux derivative lists

Version 8 retains every [v7 field](boozer-output-v7.md) and records the
meaning of the source `aphi` vector. Writers select v8 only for flux lists;
polynomial flux retains the existing versions.

`pflux_type="list"` means `aphi[i]` samples `T'(s)` at `s=i/(N-1)`.
Require at least two finite, strictly positive samples. Interpolate linearly
and integrate complete and partial trapezoids from zero to obtain `T(s)`.
Normalize the signed flux derivative with `phiedge/(2*pi*T(1))` and the
existing Jacobian sign. Sample the derivative directly on the full and half
grids; keep the half-grid lambda normalization distinct from averages of
full-grid derivatives.

The binary magic is `MCBOOZ08` and the version is 8. Append one
length-prefixed `pflux_type` string after the v7 surface-profile arrays.
The NetCDF/HDF5 schema is `magnetic-coordinate-boozer-v8`; add the string
attribute `pflux_type` alongside the existing tracing metadata. Versions
3–7 imply `power_series`, whose `aphi` coefficients define `T(s)/s`.

cuMES binary v11 carries the same selector after its input record. It always
includes a six/twelve component count after `ns, mnmax`; only twelve-family
records contain the asymmetric input extension. The reader evaluates list
pressure after reading this flux selector.

Readers that only understand v7 must reject v8 rather than interpret its
samples as polynomial coefficients. Analytic normalization tests compare
linear derivative lists with equivalent polynomials on several radial grids.
Binary tests cover both symmetries, malformed samples, and truncated records.
