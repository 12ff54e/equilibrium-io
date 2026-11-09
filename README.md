# Equilibrium data and I/O

This C++20 package reads and writes cuMES binary, HDF5, and NetCDF equilibrium
files and magnetic-coordinate Boozer exports without CUDA or either solver.
External VMEC wout files remain outside its scope. The [file compatibility contract](docs/compatibility.md) describes
the supported schemas and release rules.

```sh
git clone https://github.com/12ff54e/equilibrium-io.git
cd equilibrium-io
cmake -S . -B build \
  -DEQUILIBRIUM_USE_NETCDF=ON -DEQUILIBRIUM_USE_HDF5=ON
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/install/prefix
```

`EQUILIBRIUM_USE_NETCDF` and `EQUILIBRIUM_USE_HDF5` accept `AUTO` (the
default), `ON` (required), or `OFF`. `AUTO` skips unavailable libraries.
Binary I/O always builds. `EQUILIBRIUM_BUILD_TOOLS` and `EQUILIBRIUM_BUILD_TESTS`
default to ON for standalone builds and OFF when embedded.

HDF5 requires version 1.14.6 or newer. `AUTO` skips older versions; `ON`
rejects them. This minimum also applies to installed consumers and supplied
HDF5 targets in `equilibrium_add_io_library`.

An installed consumer selects one component:

```cmake
find_package(equilibrium_io CONFIG REQUIRED COMPONENTS io)
target_link_libraries(analysis PRIVATE equilibrium::io)
```

| Component | Interface | Dependencies |
| --- | --- | --- |
| `data` | Equilibrium and Boozer records, views, structural validation | C++20 standard library |
| `binary` | Above plus native/Boozer binary results and checkpoints | `data` |
| `io` | Above plus format selection and conversion | Enabled NetCDF/HDF5 libraries |

Selecting `data` or `binary` does not search for NetCDF, HDF5, or CUDA, even
when the installed package includes those codecs.

For a source dependency, use a pinned release or commit with FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(equilibrium_io
  GIT_REPOSITORY https://github.com/12ff54e/equilibrium-io.git
  GIT_TAG v0.2.0)
FetchContent_MakeAvailable(equilibrium_io)
target_link_libraries(analysis PRIVATE equilibrium::io)
```

The repository is public; cloning does not require a GitHub account.
`add_subdirectory(deps/equilibrium-io)` also works with a pinned submodule.
Set the codec options before adding the dependency. When embedding cuMES,
add cuMES first and link its existing `equilibrium::data`, `equilibrium::binary`,
or `equilibrium::io` targets. Do not add another copy of the package to that
build. Installed cuMES consumers reuse its `find_dependency(equilibrium_io)`.

## Boozer exports

`equilibrium/boozer.hpp` holds `BoozerResult` and its optional tracing,
magnetic-field and surface-profile groups. magnetic-coordinate computes these
records; this package validates and serializes them. It does not perform a
Boozer transform.

```cpp
#include <equilibrium/boozer_binary.hpp>
#include <equilibrium/boozer_output.hpp>

auto saved = equilibrium::read_boozer_binary("boozer.bin");
equilibrium::write_boozer_output("boozer.h5", saved.result, saved.source_path);
```

The [Boozer v8 contract](docs/boozer-output-v8.md) links the preceding schemas.
The binary component reads Boozer v7/v8 and provides `write_boozer_binary`.
The `io` component adds suffix dispatch and NetCDF/HDF5 export. `convert` and
`equilibrium-convert` recognize native equilibria and Boozer binaries. Boozer
NetCDF/HDF5 input is not supported by the C++ API. Writers emit v8 and require every record defined by that schema. They retain
source provenance and reject partial data before opening the output.

## Requested quantities

The C++ readers select physical quantities independently of file versions.
A new `ReadRequest` ignores every quantity until the caller requests it:

```cpp
#include <equilibrium/binary.hpp>

equilibrium::ReadRequest request;
using Q = equilibrium::Quantity;
using R = equilibrium::Requirement;
request[Q::GEOMETRY] = R::REQUIRED;
request[Q::CONTRAVARIANT_FIELD] = R::REQUIRED;
request[Q::GRID] = R::REQUIRED;
request[Q::PROVENANCE] = R::OPTIONAL;
auto decoded = equilibrium::read_binary("equilibrium.bin", request);
if (!decoded) throw std::runtime_error(decoded.error());
const auto& eq = decoded.value();
if (eq.has_provenance()) {
    // eq.provenance is a typed optional record.
}
const auto field = eq.contravariant_field().value();
const auto i = field.grid.index(0, 1, 0);  // radial, theta, zeta
const double b_theta = field.components[1][i];
```

`read(path, request)` and `read(path, format, request)` apply the same
contract to native binary, checkpoints, NetCDF and HDF5. Required quantities
must be present and valid. Optional quantities return usable data or a
diagnostic; ignored quantities skip allocation, decoding and validation where
the container permits. Binary readers still need framing to locate later
data. HDF5 uses the `rmncc` dataset shape for mandatory source dimensions,
even when geometry is ignored. Damage to that framing fails the read.

`DataGroup` assignments expand presets immediately; the last assignment wins.
For example, assigning `FIELDS = REQUIRED`, then `CURRENT_DENSITY = IGNORED`,
requires the magnetic Jacobian and both magnetic vectors without current
density. `GEOMETRY` as a group includes both R/Z geometry and the spectral
shift; `Quantity::GEOMETRY` includes only R/Z. The shift is native lambda or
Boozer nu. Groups have no single requirement getter after overrides; query
each quantity instead.

Validation first checks each quantity, then relations between available
quantities. Required values take precedence over conflicting optional ones.
Conflicting required quantities fail the read; conflicting optional quantities
are omitted. Bad pressure does not discard usable I/G profiles.

| Read status | Meaning |
| --- | --- |
| `AVAILABLE` | The quantity passed validation and is exposed |
| `ABSENT` | The format or file has no such quantity |
| `NOT_REQUESTED` | The caller ignored the quantity |
| `INVALID` | Decoding or validation failed; `diagnostic` explains why |
| `PARTIAL` | Group summary only: some, but not all, members are available |

Inspect `eq.read_diagnostics[Quantity::PROVENANCE]`, for example. Diagnostics
describe the read; `has_*` queries and optional members describe the object's
current contents if callers later mutate it. `has_quantity()` and named
optional accessors report individual capabilities; `has_group()` requires
every member of the preset. Coordinate dimensions remain mandatory. Native
magnetic values use the half grid; current density uses the full grid.
Grid metadata can be read without exposing or validating saved input.
Solved native profiles remain absent in the current file formats.

Boozer requests use `BoozerReadOptions::request`, which accepts quantities
and group presets. Both native and Boozer readers default to
`ReadRequest::all()`: every quantity is optional. Mark consumer prerequisites
`REQUIRED`; inspect diagnostics for optional failures. Set `options.request = {}`
to start with no requested quantities. Shared radial coordinates are always
validated. Iota is independently selectable. Magnetic data live in the single
`BoozerResult::magnetic` record and do not require tracing data.
`magnetic_data()` returns a typed optional reference.

```cpp
equilibrium::BoozerReadOptions options;
options.spectrum = equilibrium::SpectrumStorage::REAL;
options.request[equilibrium::DataGroup::GEOMETRY] = equilibrium::Requirement::REQUIRED;
options.request[equilibrium::DataGroup::TRACING] = equilibrium::Requirement::REQUIRED;
options.request[equilibrium::DataGroup::MAGNETIC] = equilibrium::Requirement::REQUIRED;
auto file = equilibrium::read_boozer_binary("boozer.bin", options);
const auto& coefficients = *file.result.real_spectrum;
```

`REAL` is the default. It preserves the requested real parity coefficients
exactly and leaves complex arrays empty. R/Z geometry and spectral nu can be
requested separately
in all three storage modes. Inactive `sin(0)` terms remain stored but contribute
nothing to complex reconstruction. `COMPLEX` reconstructs complex coefficients;
`BOTH` stores both. `has_real_spectrum()` and `has_complex_spectrum()` report
which representations are exposed. Writers use the real representation when
present, so a `REAL` read/write cycle preserves its coefficients. Writers
reject partial selections that cannot form a complete output record.

Binary reads accept unused trailing data and require valid framing for
requested quantities. `detect_binary_kind(path)`
classifies native versus Boozer input without making a capability claim.
All binary byte-order conversion belongs to this package.

Named accessors return concrete spans or views, without type erasure or copies
of numerical arrays. Field and spectral views provide their storage indices;
`full_s()` and `half_s()` provide source radial coordinates. `axis_geometry()`
evaluates R/Z at a physical toroidal angle. `b2_jacobian()` exposes the legacy
`b2j00/(4*pi*pi)` profile without allocating a normalized array; it is distinct
from the signed `magnetic_jacobian()` used with magnetic components.
Native `axis_geometry()` derives its view from R/Z and GRID; the selectable
`AXIS_GEOMETRY` quantity refers to the separate stored Boozer axis record.

Keep the owner alive and its buffers unchanged while using a view. Borrowing
from temporary owners is disabled. For caller-built data, use
`validate_quantity(as_view(owner), quantity)` for each needed quantity and
its shared metadata before taking views. Native geometry indexing also
requires a valid `GRID`; `validate()` checks a complete native record. Native
validation is header-only in `equilibrium/validation.hpp`; Boozer validation
and cross-quantity checks use `equilibrium/boozer_validation.hpp` and link
`equilibrium::binary` (or `equilibrium::io`). Consumers still own Fourier
reconstruction, interpolation and application-specific numerical checks.

## Standalone WebAssembly converter

Use Emscripten 6.0.9 and CMake 3.31.6, matching the tested CI toolchain.
Activate the SDK and put that CMake version on `PATH` before configuring.

```sh
emcmake cmake -S . -B build-wasm -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-wasm
ctest --test-dir build-wasm --output-on-failure
```

This builds `equilibrium-formats.wasm`, the codec libraries and tests, and
`licenses/` notices. CMake downloads pinned zlib, HDF5 and NetCDF sources and
builds them for Wasm. It does not search for native scientific libraries.
Both codec options must be enabled to build the full converter. Set both to
`OFF` for binary-only library tests, or `EQUILIBRIUM_BUILD_TOOLS=OFF` to omit
the converter. The usual `equilibrium::data`, `binary`, and `io` targets remain
available to source consumers.

```js
import {createFormatsCpu} from './js/formats_cpu.mjs';
const formats = await createFormatsCpu(wasmBytes);
const hdf5Bytes = formats.convert(equilibriumBytes, 'h5');
```

Ship `cpu_wasm.mjs` and `cpu_files.mjs` beside `formats_cpu.mjs`, plus the
license notices. Conversion uses unshared Wasm memory and copies owned bytes
through WasmFS. The same module handles native and Boozer exports.

Embedded Emscripten builds can include `cmake/EquilibriumWasm.cmake` and call
`equilibrium_add_wasm_converter(target TESTS)`. Omit `TESTS` to omit the
converter fixtures. The target's `EQUILIBRIUM_LICENSE_FILES` property lists the
notices to stage with the Wasm artifact. The helper owns all source and codec
paths; it needs no cuMES or magnetic-coordinate checkout.

## Python readers

Install from this checkout with `python3 -m pip install .`. Binary and checkpoint
reads need only NumPy. Use `pip install '.[netcdf]'`, `'.[hdf5]'`, or `'.[all]'`
for SciPy NetCDF and h5py HDF5 support. No cuMES or plotting package is required.

```python
from equilibrium_io import load_state
ns, mnmax, families, input_params, name = load_state("equilibrium.bin")
```

`load_state` preserves the plotting reader's tuple and mode-major NumPy arrays.
It requires saved input metadata. Format checks raise `ValueError` rather than
terminating the process; ordinary file and truncated-struct errors remain
catchable Python exceptions. `load_state` reads native cuMES formats only.

The Python API also reads the geometry and field-strength subset of Boozer
v7/v8 binary, NetCDF and HDF5 files:

```python
from equilibrium_io import load_boozer
data = load_boozer("boozer.bin")
print(data.s, data.iota, data.b.shape)
```

`BoozerData` contains real Fourier geometry families and `b` in
`[surface, zeta, theta]` order. Its `alpha` property is the source field-period
angle; `theta_b` is the Boozer poloidal angle. This analysis API does not retain
all optional records and cannot rewrite a complete file without data loss.
Malformed files raise `ValueError`; filesystem errors raise `OSError`.
Rendering and coordinate-grid construction belong to the plotting consumer.
The Python reader needs the same optional codec extras listed above.

## Spectral payload inspection

`equilibrium/spectral_payload.hpp` exposes `read_spectral_payload(path)` and
`read_spectral_state(path)` through `equilibrium::binary`. These inspect only
the spectral portion of result binaries v10/v11. They accept missing or
unsupported trailing metadata and do not require a derived-field block.
`read_spectral_payload` returns the exact stored bytes, including signed zero
and NaN payload bits; it does not decode and rewrite them. Unversioned and older files are rejected. The per-family resource cap applies. Use `read` or `read_binary` for complete equilibria and
checkpoints.

For browser diagnostics, `js/binary_payload.mjs` exports
`scientificPayloadRanges(bytes)` and `digestScientificPayload(bytes)` for
result binaries v10/v11 with derived fields. The latter uses Web Crypto SHA-256
and excludes headers and provenance from its state and field digests. Import
this module directly; it has no npm dependencies.

## Alternate CMake toolchains

Source consumers that need several compiler variants can include
`cmake/EquilibriumTargets.cmake`. `equilibrium_add_binary_library(target)` builds
the binary codecs; `equilibrium_add_io_library(target binary_target
NETCDF_TARGET netcdf_target HDF5_TARGET hdf5_target)` adds the selected scientific
codecs. Omit either dependency argument to disable that codec. These helpers
use supplied targets without searching the host system. Consumers set thread,
exception, and sanitizer flags on the resulting targets. Ordinary consumers
should use the package's `equilibrium::binary` and `equilibrium::io` targets.

## Development and releases

Run the tests with both codecs enabled, each enabled alone, and both disabled.
The package tests also install and compile external data, binary, and full-I/O
consumers. CI runs this matrix without CUDA, plus ASan/UBSan, binary-only WebAssembly tests and a full Wasm converter check.

The suite covers historical binary layouts, checkpoints, asymmetric input,
malformed container shapes and restart offsets, and atomic publication failures.
cuMES keeps its configuration adapter, device transfer, derived-field physics,
solver/CLI, plot rendering, and browser integration tests. Python reader tests
and JavaScript payload tests live here. Set `EQUILIBRIUM_PYTHON_TESTS=ON` or
`EQUILIBRIUM_JS_TESTS=ON` to include them in CTest; CI also tests installed Python
packages outside the source tree, with NumPy alone and with all codec extras.

Enable `EQUILIBRIUM_SANITIZER_TESTS=ON` to build ASan/UBSan copies of the codecs
and their tests. These copies are not installed. CTest writes fixtures beneath
the build directory. A binary-only Emscripten build runs the same codec tests
through Node; installed-consumer checks and Linux filesystem fault probes run
only in native builds.

Update the CMake version and changelog before tagging a release. Treat C++ API
compatibility and file compatibility separately: an API change does not permit
silently changing a persisted schema. Follow the
[compatibility contract](docs/compatibility.md), commit and publish the tested
library revision, then update consumer submodule pins. Consumers should cleanly
rebuild when public C++ layouts change.

The initial library history was split from cuMES commit
`35fc387` (`lib/equilibrium`). Earlier codec history remains in cuMES.

## Read, inspect, and write

```cpp
#include <equilibrium/io.hpp>
#include <stdexcept>

auto result = equilibrium::read("equilibrium.h5");
if (!result) throw std::runtime_error(result.error());
const auto& eq = result.value();
const auto view = equilibrium::as_view(eq);
const auto r = view.families[equilibrium::EquilibriumSnapshot::RMNCC];
// r[surface + mode * view.ns] borrows the owner's coefficients.
if (view.fields) {
    const auto sqrtg = view.fields->half[equilibrium::EquilibriumSnapshot::SQRTG];
    // Use sqrtg while eq and its buffers remain alive.
}
auto status = equilibrium::write("copy.nc", view);
if (!status) throw std::runtime_error(status.error());
```

`Equilibrium` owns numerical vectors. `EquilibriumView` instantiates the same
`BasicEquilibrium` template with spans and metadata references. Copying or
moving an owner is safe; obtain a fresh view afterward. A view does not extend
the owner's lifetime, and resizing an owner's vectors can invalidate it.
`as_view` rejects temporary owners.

`fields`, `profiles`, `grid`, `input`, and `provenance` are optional groups.
Absence differs from a present array of zeros. `validate(view)` checks shape
and metadata consistency; it does not establish physical validity or solver
convergence. Partial field or profile groups fail validation.

Current files do not persist solved radial profiles. Reads leave `profiles`
absent. Writes reject populated profiles unless the caller supplies
`WriteOptions{ProfileWritePolicy::OMIT}`. Input profile prescriptions in
`input` are distinct from these solved profiles.

Native binary readers accept result v10/v11 and checkpoint v8/v9.
Writers emit result v11 and checkpoint v9 for both symmetries and flux kinds.
No reader falls back to an older format or weaker request after an error.
Unused trailing bytes are allowed. To read only spectral geometry, request
`DataGroup::GEOMETRY` and leave other quantities ignored.

Checkpoints provide spectral state and saved input, without scientific fields
or run outcomes. Read them with `read_binary` or `read`; converting one to a
result requires the caller to supply run provenance. `Reader::read` returns
an `Equilibrium` and accepts the same request; it has no report out-parameter.

`source` records the original container, version, and checkpoint status.
HDF5/NetCDF use logical schema version 1, independent of binary versions.

Input detection checks file signatures. NetCDF-4 shares HDF5's signature:
use a `.nc` suffix or `read(path, OutputFormat::NETCDF)` for NetCDF-4 files.
Output suffixes select `.bin`, `.nc`, `.h5`, or `.hdf5`. Unsupported formats
and missing codecs return errors. Writers publish through a temporary file
in the destination directory. Failures before rename preserve the destination.
A directory-sync error after rename means replacement occurred but durability
could not be confirmed.

HDF5 writers build a file image in memory, close HDF5, then write the image
through the checked publication path. This keeps storage-write failures out
of HDF5 cleanup, which can leave the library unusable even in newer versions.
Image extraction briefly holds two serialized images plus HDF5 caches; one
image remains during the disk write. Allow for this extra memory when writing
large files. This protects against storage failures, not allocation failure.

```sh
equilibrium-convert input.bin output.h5
```

Conversion preserves scientific values and provenance. It does not solve,
interpolate, normalize, or fill missing data.

## Compose a view from existing owners

An analysis program can fill `EquilibriumView` directly from its arrays:

```cpp
equilibrium::EquilibriumView view;
view.ns = ns;
view.mnmax = mpol * (ntor + 1);
view.families = {rcc, zsc, lsc, rss, zcs, lcs}; // spans into separate owners
view.grid = equilibrium::GridMetadata{mpol, ntor, nfp, ntheta, nzeta};
auto status = equilibrium::validate(view);
```

This needs no cuMES configuration. Attach optional groups only when they are
available. `as_equilibrium_view(snapshot)` also adapts the legacy host
snapshot. The cuMES header `cumes/io/equilibrium_adapter.hpp` composes a
snapshot, saved input, optional solved profiles, and optional run report.
Device capture, paired-float conversion, radius-reference restoration, and
derived-field computation remain in cuMES.

meow may borrow this data while retaining its owning solver outcome, cache,
convergence checks, and terminal-correction evidence. Tangents and cotangents
have different meanings from physical equilibria and must not pass through
ordinary result publication.

## Numerical conventions

Spectral arrays store physical doubles, with `surface + mode * ns` indexing
and `mode = m * (ntor + 1) + n` for folded `n >= 0`. Symmetric families are
`Rcc, Zsc, Lsc, Rss, Zcs, Lcs`; asymmetric states append
`Rsc, Zcc, Lcc, Rcs, Zss, Lss`. The letters name cosine/sine factors in poloidal
and toroidal angles. Physical toroidal mode number is `n*nfp`. Theta covers
`[0, 2*pi)`; physical zeta covers one field period `[0, 2*pi/nfp)`.
R and Z use metres; lambda is an angular displacement in radians.

The full radial grid uses `s=j/(ns-1)`; the half grid uses
`s=(j+1/2)/(ns-1)`. Real-space indexing is
`theta + ntheta * (zeta + nzeta * surface)`. Half-grid fields are
`sqrtg, bsups, bsupu, bsupv, bsubs, bsubu, bsubv`; full-grid fields are
`jsups, jsupu, jsupv, jsubs, jsubu, jsubv`. Superscript components are
contravariant; subscript components are covariant in `(s, theta, zeta)`.
These are coordinate components, not Cartesian components: with dimensionless
coordinates, magnetic units are T/m and T*m respectively; current units are
A/m³ and A/m. The signed coordinate Jacobian has units m³.

Solved profiles use the half grid. Flux derivatives are Wb per unit `s` in
the VMEC wout convention. Iota is dimensionless; `I=<B_theta>` and
`G=<B_zeta>` use T*m and normalized angular quadrature. No adapter changes
these conventions.

## cuMES compatibility

The original cuMES headers, factory calls, abstract `Writer` interface,
checkpoint calls, and solver APIs remain available. They forward to one codec
implementation and borrow numerical arrays when writing. Existing consumers
must rebuild: this extraction changes C++ type identities and is not an ABI
compatibility promise. `RunReport` now extends `RunProvenance`; direct field
assignment still works, but designated aggregate initialization of inherited
members must be replaced by assignment or base initialization.
