# Equilibrium data and I/O

This C++20 package reads and writes cuMES binary, HDF5, and NetCDF equilibrium
files without CUDA or the cuMES solver. Its source directory can be copied and
built on its own. VMEC wout and Boozer files are outside its scope.

```sh
cmake -S lib/equilibrium -B build-equilibrium \
  -DEQUILIBRIUM_USE_NETCDF=ON -DEQUILIBRIUM_USE_HDF5=ON
cmake --build build-equilibrium
ctest --test-dir build-equilibrium --output-on-failure
cmake --install build-equilibrium --prefix /your/install/prefix
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
