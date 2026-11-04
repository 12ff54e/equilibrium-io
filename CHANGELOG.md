# Changelog

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
