# Equilibrium data and I/O

This C++20 package reads and writes cuMES binary, HDF5, and NetCDF equilibrium
files without CUDA or the cuMES solver. VMEC wout and Boozer files are outside
its scope. The [file compatibility contract](docs/compatibility.md) describes
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

Disable either codec when its development libraries are unavailable. Binary
I/O always builds. `EQUILIBRIUM_BUILD_TOOLS` and `EQUILIBRIUM_BUILD_TESTS`
default to ON for standalone builds and OFF when embedded.

An installed consumer selects one component:

```cmake
find_package(equilibrium_io CONFIG REQUIRED COMPONENTS io)
target_link_libraries(analysis PRIVATE equilibrium::io)
```

| Component | Interface | Dependencies |
| --- | --- | --- |
| `data` | Records, views, structural validation | C++20 standard library |
| `binary` | Above plus binary results and checkpoints | `data` |
| `io` | Above plus format selection and conversion | Enabled NetCDF/HDF5 libraries |

Selecting `data` or `binary` does not search for NetCDF, HDF5, or CUDA, even
when the installed package includes those codecs.

For a source dependency, use a pinned release or commit with FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(equilibrium_io
  GIT_REPOSITORY https://github.com/12ff54e/equilibrium-io.git
  GIT_TAG v0.1.0)
FetchContent_MakeAvailable(equilibrium_io)
target_link_libraries(analysis PRIVATE equilibrium::io)
```

The repository is public; cloning does not require a GitHub account.
`add_subdirectory(deps/equilibrium-io)` also works with a pinned submodule.
Set the codec options before adding the dependency. When embedding cuMES,
add cuMES first and link its existing `equilibrium::data`, `equilibrium::binary`,
or `equilibrium::io` targets. Do not add another copy of the package to that
build. Installed cuMES consumers reuse its `find_dependency(equilibrium_io)`.

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
catchable Python exceptions. The package reads cuMES formats only.

## Spectral payload inspection

`equilibrium/spectral_payload.hpp` exposes `read_spectral_payload(path)` and
`read_spectral_state(path)` through `equilibrium::binary`. These inspect only
the spectral portion of result binaries v1-v11. They accept missing or
unsupported trailing metadata and do not require a derived-field block.
`read_spectral_payload` returns the exact stored bytes, including signed zero
and NaN payload bits; it does not decode and rewrite them. Legacy unversioned
files require `allow_legacy=true` and an exact file size. The usual per-family
resource cap applies. Use `read` or `read_binary` for complete equilibria and
checkpoints.

For browser diagnostics, `js/binary_payload.mjs` exports
`scientificPayloadRanges(bytes)` and `digestScientificPayload(bytes)` for
result binaries v8-v11 with derived fields. The latter uses Web Crypto SHA-256
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
consumers. CI runs this matrix without CUDA, plus ASan/UBSan and binary-only
WebAssembly checks.

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

Complete reads are the default. They require saved input and, for result
files, valid run provenance. No reader retries in a weaker mode after an
error. Use `read(path, ReadMode::STATE_ONLY)` explicitly for historical state
recovery: it skips the provenance trailer and leaves input, grid metadata,
and provenance absent. Binary versions 1–11 and checkpoint versions 1–9
support state recovery. Complete reads reject the obsolete JSON input in
binary v9 and checkpoint v7, and records without complete saved input.
Complete binary reads reject invalid tags, truncation, and trailing bytes.

Checkpoints return spectral state and, in complete mode, saved input. They
contain neither scientific fields nor run outcomes. Converting a checkpoint
to a result therefore requires the caller to supply run provenance.
`source` records the original container, version, and checkpoint status.
HDF5/NetCDF use logical schema version 1, independent of binary versions.

Input detection checks file signatures. NetCDF-4 shares HDF5's signature:
use a `.nc` suffix or `read(path, OutputFormat::NETCDF)` for NetCDF-4 files.
Output suffixes select `.bin`, `.nc`, `.h5`, or `.hdf5`. Unsupported formats
and missing codecs return errors. Writers publish through a temporary file
in the destination directory. Failures before rename preserve the destination.
A directory-sync error after rename means replacement occurred but durability
could not be confirmed.

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
