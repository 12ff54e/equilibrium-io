# File compatibility contract

The package owns cuMES result/checkpoint codecs and magnetic-coordinate Boozer
codecs. It does not depend on the cuMES solver.

| Container | Accepted versions | Written version |
| --- | --- | --- |
| Binary result, `CUMES001` | 10, 11 | 11 |
| Checkpoint, `CUMECKP1` | 8, 9 | 9 |
| NetCDF / HDF5 native result | Logical schema 1 | 1 |
| Boozer binary, `MCBOOZ` | 7, 8 | 8 |
| Boozer NetCDF / HDF5 | Python analysis reads 7, 8; C++ export only | 8 |

Only the latest two versions of each versioned format are supported. Older,
unversioned, and future versions fail before payload decoding. This applies
to full reads, selected reads, Python analysis, and C++/JavaScript payload
inspection. Native NetCDF/HDF5 still has one logical schema. NetCDF-4 needs a
`.nc` suffix or explicit codec selection because it shares HDF5's signature.

All current writers emit the latest schema, regardless of symmetry or flux
kind. Boozer v8 always includes its symmetry flag, tracing, magnetic fields,
angle maps, surface profiles, and trailing `pflux_type`. The selector accepts
`power_series` and `list`; v7 implies `power_series`. A writer rejects an
incomplete result before opening its output. `source_format_version` is source
provenance, not the Boozer container version.

Result v10 and checkpoint v8 always include an asymmetric input extension,
even for six-family state. Result v11 and checkpoint v9 include that extension
only for twelve-family state, then append the flux selector.

## Requests and ownership

All C++ reads use `ReadRequest`. No-options reads use `ReadRequest::all()`,
which treats every quantity as optional. Consumers mark their prerequisites
`REQUIRED`. Optional quantities fail independently and retain diagnostics;
group assignments expand presets, and later assignments win. Coordinate
framing remains mandatory. Unused trailing bytes are allowed.

Group capability queries require every member; group diagnostics can report
`PARTIAL`. Writers reject selections missing required output data. Conversion
requests all quantities and rejects invalid data before writing. It preserves
values and provenance without solving or resampling. Native files contain no
solved radial profiles; input prescriptions do not replace them.

Boozer reads default to `SpectrumStorage::REAL`, preserving stored real
coefficients, including inactive terms. `COMPLEX` and `BOTH` remain runtime
options. Magnetic arrays have one owner in `BoozerResult::magnetic`.

The C++ API has no backward compatibility layer. `ReadMode`, `read_checkpoint`,
the `Reader` report out-parameter, `allow_legacy`, and nested magnetic storage
have been removed. Consumers must migrate and rebuild. Producer snapshot/view
adapters remain available for writers.

## Data and publication

Preserve binary64 storage, family order, radial staggering, layouts, units,
raw harmonic order, duplicates, explicit zeros, folded boundaries, profile
selectors, and provenance. Scientific convention labels such as
`mixed-grid-v2` and `real-parity-v3` describe mathematics, not admission versions.
The linked older Boozer documents define blocks retained in current layouts.

Python readers raise `ValueError` for format errors and retain their analysis
array layouts. Spectral-payload inspection reads only the result's spectral
bytes, without requiring later fields or metadata or checking finiteness.
Comparison tools own numerical acceptance and bitwise-equality policy.

Writers rename a completed temporary file in the destination directory.
Pre-rename failures preserve the destination. A directory-sync failure after
rename reports unconfirmed durability after replacement.
