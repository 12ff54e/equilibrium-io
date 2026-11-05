# Asymmetric Boozer output format version 4

Version 4 extends the [version-3 contract](boozer-output-v3.md) to inputs
without stellarator symmetry. These are legacy formats without tracing data;
new transforms use [version 7](boozer-output-v7.md) for both symmetries.
The binary reader accepts versions 3 through 7. The transform accepts cuMES
binary version 8
(six families) and version 10 (twelve families plus typed asymmetric input).
The obsolete JSON-based cuMES version 9 is not supported.

The owning and borrowed equilibrium types default to six spectral families.
Resize `families` to `ASYMMETRIC_SPECTRAL_FAMILY_COUNT` for asymmetric data;
`lasym()` reports this extent. Input family order is
`Rcc, Zsc, Lsc, Rss, Zcs, Lcs, Rsc, Zcc, Lcc, Rcs, Zss, Lss`, with
`surface + mode*ns` inside each family. All twelve families are required.
`BoozerResult::lasym` selects the output representation.

## Coordinates and amplitudes

The angular map, units, radial interpolation, axis exclusion and Fourier
normalization are unchanged. The output is uniform in `theta_b` and unchanged
`alpha=nfp*zeta`; `alpha_b=alpha+nfp*nu`, where `nu` is in physical toroidal
radians. It is a mixed grid, not a grid uniform in both Boozer angles.

Each of `R`, `Z` and `nu` has all four real reconstruction amplitudes:

```text
f = sum [fcc cos(m theta_b) cos(n alpha)
       + fss sin(m theta_b) sin(n alpha)
       + fsc sin(m theta_b) cos(n alpha)
       + fcs cos(m theta_b) sin(n alpha)].
```

The field prefixes are `r`, `z`, and `nu`; for example, `numncc` stores the
cosine-cosine part of `nu`. Modes remain nonnegative, ordered by `n` then `m`.
Arrays use `[surface][mode]`. No additional FFT normalization applies when
reconstructing a field. Sine factors vanish when their corresponding mode is
zero. The constant `Z` coefficient is retained, including vertical offsets.

The internal complex spectrum must remain Hermitian for each real field.
Serialization verifies reconstruction against that spectrum using the same
`1e-10 * max(1, abs(coefficient))` bound as version 3. It does not impose even
parity on `R` or odd parity on `Z`/`nu` for asymmetric inputs. No symmetry
averaging is performed. Signed Jacobian, positive field strength, monotone
angular maps, nonzero flux/zero modes and resonance checks remain active.

## Binary layout

Integers are little-endian signed int32; real values are little-endian
IEEE-754 binary64. Strings are an int32 byte length followed by UTF-8 bytes.

1. Eight bytes `MCBOOZ04`, followed by int32 version `4`. This version
   identifies the full asymmetric representation; there is no extra flag.
2. Three strings: coordinate convention, Fourier convention, source path.
   The coordinate string is unchanged from version 3. The Fourier string
   begins `real-parity-v4` and describes all four product terms.
3. Thirteen int32 values: `source_format_version`, `source_ns`,
   `source_ntheta`, `source_nzeta`, `source_mpol`, `source_ntor`, `nfp`,
   `first_surface`, `ntheta`, `nzeta`, `mmax`, `nmax`, `radial_order`.
4. One binary64 `resonance_tolerance`.
5. Binary64 arrays `s`, `iota`, each with `source_ns-first_surface` entries.
6. Int32 arrays `mode_m`, `mode_n`, each with `(mmax+1)*(nmax+1)` entries.
7. Twelve binary64 coefficient arrays, each of size `surfaces*modes`, in
   order: `rmncc`, `rmnss`, `zmnsc`, `zmncs`, `numnsc`, `numncs`,
   `rmnsc`, `rmncs`, `zmncc`, `zmnss`, `numncc`, `numnss`.
8. Binary64 grids `B`, `sqrt_g_b`, each ordered
   `[surface][alpha][theta_b]`, followed by `b2j00[surface]`.

Array sizes are determined by the header; arrays have no length prefixes.
The invariant `B^2 * sqrt_g_b = b2j00 / (4*pi^2)` preserves the signed source
Jacobian orientation, as in version 3.

## NetCDF and HDF5

Both containers use schema `magnetic-coordinate-boozer-v4` and the root/global
integer attribute `lasym=1`. All twelve coefficient variables/datasets use
the names and extents above. Version-3 dimensions, other metadata and physical
fields are unchanged. Coefficients remain ordinary real rank-2 double arrays.

## Validation

Manufactured tests cover each additional input family, physical toroidal
derivatives with `nfp>1`, zero-mode limits, asymmetric Boozer shift harmonics,
and complete axisymmetric/3-D transforms with an analytic identity map. Output
tests reconstruct analytic real fields from the binary spectrum, check each
new NetCDF/HDF5 family and reject symmetry loss, non-Hermitian spectra,
incomplete inputs and truncated records. CUDA geometry and Fourier operations
are compared with independent host evaluation. The same tests run under Wasm.
