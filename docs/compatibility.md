# File compatibility contract

The package owns cuMES result/checkpoint codecs and magnetic-coordinate Boozer
result codecs. It does not depend on the
cuMES solver and does not interpret another program's format as cuMES data.

| Container | State recovery | Complete record |
| --- | --- | --- |
| Binary result, `CUMES001` | Versions 1–11 | Typed saved input in v3–8 and v10–11; full provenance and structural validation required |
| Checkpoint, `CUMECKP1` | Versions 1–9 | Typed saved input in v2–6 and v8–9; no run provenance |
| NetCDF / HDF5 result | Logical schema 1 | Required scientific structure, saved input, and run provenance |
| Boozer binary, `MCBOOZ` | Versions 3–8 | Optional groups follow the stored version; source path retained |
| Boozer NetCDF / HDF5 | Export only | Existing v3–v8 schemas and real-parity normalization |

Binary v9 and checkpoint v7 used an obsolete JSON input extension. Explicit
state-only reads recover their spectral state; complete reads reject their
input. Older records that lack the input required for analysis also need an
explicit state-only read. A failed complete read never falls back to that mode.

Current writers retain the existing version selection: binary v8 and
checkpoint v6 for symmetric polynomial-flux input, binary v10 and checkpoint
v8 for asymmetric polynomial-flux input, and binary v11 and checkpoint v9 for
list flux input. NetCDF output uses the existing classic schema; readers also
accept its NetCDF-4 representation. NetCDF-4 needs a `.nc` suffix or explicit
codec selection because it shares HDF5's signature.

Spectral state remains double on disk. Preserve family order, radial
staggering, mode and point layouts, units, raw harmonic ordering, duplicates,
explicit zeros, folded boundaries, profile selectors, and provenance. Container
conversion preserves decoded values and metadata without solving or resampling.
The README defines the scientific array conventions.

Missing groups remain absent. Current result files can contain the complete
scientific-field block; checkpoints cannot. Neither format stores solved
radial profiles. Writers reject populated profiles unless the caller explicitly
chooses omission. Input profile prescriptions do not replace solved profiles.

A format extension must introduce a version or schema rule that old readers
can reject safely. It must not reinterpret bytes under an existing version or
invent missing physical data. Keep historical state recovery separate from
complete-record admission. Add regression fixtures for new formats and retain
coverage for the supported historical versions. Report intentional reader
compatibility changes in the changelog.

Writers publish by renaming a completed temporary file in the destination
directory. Pre-rename failures preserve the destination. A directory-sync
failure after rename reports unconfirmed durability after replacement.

The initial extraction passed byte comparisons for five binary fixtures and
full-precision schema/value comparisons for ten NetCDF/HDF5 fixtures, covering
symmetry, asymmetry, free-boundary provenance, and list profiles. This library
now owns the historical-reader, malformed-container, codec, and installed-package
regression suites. cuMES retains adapter, derived-field, solver, and browser
integration checks.

Requested-group reads provide a separate admission contract: required,
optional and ignored groups have independent validation and diagnostics.
See [Requested data groups](../README.md#requested-data-groups). These reads
use capabilities rather than consumer version gates and permit unused
trailing binary data. Existing no-options APIs retain their strict behavior.

## Analysis readers

The Python `equilibrium_io.load_state` API retains the former cuMES plotting
reader's tuple, family ordering, NumPy array layout, and input dictionaries.
It does not import plotting modules. Format rejection raises `ValueError`;
callers can catch errors instead of handling process termination.

C++ spectral-payload inspection accepts a result header and its complete state
payload without reading later fields or provenance. This is a distinct contract
from `ReadMode::STATE_ONLY`, which also reads any stored derived-field block.
Byte inspection preserves every stored bit and does not check finiteness.
Legacy files remain opt-in. Comparison tools retain responsibility for numerical
acceptance and bitwise-equality policy.

Boozer data types moved from magnetic-coordinate to this package in 0.2.
magnetic-coordinate retains source aliases, but consumers must rebuild because
the C++ ABI changed. Flux lists select Boozer v8; surface profiles select v7;
older results keep their older schema. Default complex reads may round spectral
parity pairs through the complex/real normalization. Select
`SpectrumStorage::REAL` or `BOTH` to retain the stored real coefficients exactly.
Writers prefer retained real arrays. Sampled fields, profiles, angle maps and
metadata remain exact. No mathematical transform or resampling occurs during
container conversion.
