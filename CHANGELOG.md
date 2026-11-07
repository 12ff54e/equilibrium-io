# Changelog

## Unreleased

- Select individual quantities with group assignment presets. Validate optional
  values independently, including pressure versus I/G and R/Z versus spectral
  shifts. Report partial groups without hiding usable members.
- Add named borrowed views, shared native/Boozer validation, field and spectral
  indexing, radial/axis helpers, and normalized legacy Boozer B²J access.
  These APIs use concrete types and do not copy numerical arrays.
- Add required, optional and ignored group requests for native binary,
  checkpoint, NetCDF, HDF5 and Boozer binary reads, with per-group diagnostics
  and instance capability queries. Existing no-options reads remain strict.
- Retain original real Boozer coefficients with runtime REAL/COMPLEX/BOTH
  selection. REAL reads avoid complex arrays; writers prefer retained real
  values. Optional profiles no longer constrain independent selected groups.
- Handle native binary scalar and array byte order inside the codecs.
- Expose binary kind detection for consumers without file-version gates.
- Extend public data records and the Reader interface; C++ consumers must
  rebuild against these headers.

- Require HDF5 1.14.6 or newer, including supplied codec targets and installed
  consumers. Native `AUTO` builds skip older HDF5 libraries.
- Build HDF5 output in memory before publishing it, so storage-write failures
  preserve the destination, remove temporary files, and leave HDF5 usable.
  Image extraction briefly requires two serialized images plus HDF5 caches.
- Check HDF5 dataset and attribute closes and release handles on early exits.
- Test repeated storage failures, resource cleanup, library shutdown, and a
  later successful write for both state and Boozer output.

## 0.3.0 — 2026-11-05

- Move the Python Boozer v3–v8 binary, NetCDF and HDF5 analysis readers from
  cuMES into this package. `load_boozer` returns geometry, iota and field
  strength; it does not retain every optional record for lossless rewriting.
- Raise ordinary Python exceptions in the library. Plotting applications
  handle command-line exits.

## 0.2.0 — 2026-11-05

- Own Boozer data records, v3–v8 binary reads, and binary/NetCDF/HDF5 exports.
  Consumers of the moved C++ types must rebuild. Boozer scientific-container
  reads remain unsupported.
- Publish Boozer files through the same checked temporary-file/rename path as
  native equilibria, preserving destinations on failures before publication.
- Own the standalone Wasm converter, pinned scientific dependencies, filesystem
  bridge, JavaScript loader, notices, and conversion tests.
- Accept AUTO/ON/OFF codec policies; AUTO skips missing native libraries.

- Add an installable Python reader and its cross-language fixture tests.
- Add spectral-only C++ reads that preserve exact payload bytes, and a
  JavaScript module for scientific-payload ranges and digests.
- Share codec target construction across native, sanitizer, and Wasm builds.

- Move the cuMES codec regression tests into this repository: historical binary
  layouts, checkpoints, asymmetric provenance, malformed containers, and atomic
  publication failures. They run without cuMES or CUDA.
- Add ASan/UBSan codec tests and binary-only WebAssembly CI.
- Make the repository public.

## 0.1.0 — 2026-11-05

- Extract the cuMES equilibrium records and codecs into an independent C++20
  package with data, binary, and optional HDF5/NetCDF components.
- Provide owning records, borrowed views, explicit missing groups, complete
  and state-only reads, atomic result publication, and a converter.
- Preserve existing cuMES file schemas. Solved radial profiles remain
  in-memory data and require explicit omission when writing current formats.
- Test standalone builds and installed consumers without CUDA.

Existing cuMES consumers must rebuild against the new C++ types. See the
README for source compatibility limits and the cuMES adapter.
