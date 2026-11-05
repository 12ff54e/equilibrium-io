# Boozer output format version 5

Version 5 introduced signed flux and axis data for both symmetries.
New transforms write [version 7](boozer-output-v7.md), which adds magnetic
components in the Boozer basis, untruncated coordinate maps, pressure,
and the covariant coefficients `I` and `G`. The non-axis
mixed grid, real Fourier amplitudes,
signed Jacobian, and angular conventions of [v3](boozer-output-v3.md) and
[v4](boozer-output-v4.md) are unchanged. Version 5 adds the signed flux
profiles and magnetic-axis data needed for downstream field reconstruction.
The transform does not perform particle tracing or choose the consumer's
interpolation or axis regularization.

`BoozerResult::tracing` is populated by every new transform. The binary reader
also accepts v3/v4, leaving `tracing` absent: their `|B|`, iota and Jacobian
cannot distinguish a field from its reversal. Rewriting a legacy result
without tracing data preserves its legacy format; no sign or axis is invented.

## Signed flux and radial coordinate

`s` is the **source radial coordinate** `j/(source_ns-1)`. It is proportional
to toroidal flux only when the input `aphi` profile is constant. Version 5
retains both the original flux inputs and the profiles recovered from the
converged field:

| Name | Extent | Meaning |
| --- | --- | --- |
| `phiedge` | scalar | Original signed cuMES edge toroidal flux, Wb |
| `aphi` | `aphi_count` | Original dimensionless power-series coefficients |
| `toroidal_flux_derivative` | `source_ns` | Signed `phip = <sqrt(g) B^zeta>` on all full-grid surfaces |
| `poloidal_flux_derivative` | `source_ns` | Signed `chip = <sqrt(g) B^theta>` on all full-grid surfaces |

Here `zeta` is the **physical** source toroidal angle, `theta` the native
poloidal angle, and the brackets are uniform angular averages on each native
half surface. The same two- or four-point radial interpolation used by the
transform places these averages at `s=j/(source_ns-1)`, including `j=0`.
The endpoints are explicit polynomial extrapolations, not extra solver field
evaluations. The units are Wb per unit dimensionless `s`, for flux divided by
`2*pi`. No extra `nfp` factor or absolute value applies. The exported non-axis
`iota` equals `chip/phip`; the profiles also supply the axis value of this ratio.

For the original cuMES input convention, define

```text
T(s) = sum_k aphi[k] * s^(k+1),
psi_input(s) = -phiedge * T(s) / (2*pi*T(1)).
```

The minus sign is cuMES's flux orientation convention. `psi_input(0)=0`;
its derivative is the prescribed toroidal flux derivative. Use the exported
field-derived `phip` with the exported Jacobian to reconstruct the converged
field. Both inputs and recovered profiles are retained so a consumer need
not assume they agree exactly at finite numerical resolution.

For example, away from the axis, write the embedding in physical space as
`x(s,u,a)=(R*cos(a/nfp), R*sin(a/nfp), Z)`, with `u=theta_b`, `a=alpha`.
The stored `nu` obeys `zeta_b=a/nfp+nu(s,u,a)`. The usual straight-field-line
contravariant representation gives

```text
B = (phip / sqrt_g_b) *
    [iota*x_u + (1 - iota*nu_u)/(1/nfp + nu_a) * x_a].
```

Subscripts denote partial derivatives on the mixed grid. The signed
`sqrt_g_b` is the Jacobian for **physical** Boozer toroidal angle, despite
being sampled at uniform source `alpha`. Thus this expression includes the
physical-angle factor once. Geometry, shift, signed flux, and the sampled
`|B|` and Jacobian provide the inputs for a consumer's consistent field and
derivative interpolation. Truncation/interpolation errors remain the
consumer's responsibility. Reversing all covariant and contravariant source
field components negates `phip`, `chip`, and `phiedge`; it leaves `|B|`,
iota, geometry, shift and Jacobian unchanged.

## Magnetic axis and its field limit

`axis_r_cos`, `axis_r_sin`, `axis_z_cos`, `axis_z_sin` each contain
`source_ntor+1` real reconstruction amplitudes in metres, ordered `n=0..source_ntor`:

```text
R_axis(alpha) = sum_n [axis_r_cos[n]*cos(n*alpha)
                    + axis_r_sin[n]*sin(n*alpha)],
Z_axis(alpha) = sum_n [axis_z_cos[n]*cos(n*alpha)
                    + axis_z_sin[n]*sin(n*alpha)].
```

These are the converged source `m=0` coefficients, including asymmetric
radial sine and vertical cosine terms. They are independent of output
angular resampling and spectral truncation. They are not the input axis
guess, and do not include the nonphysical `m>0` entries that cuMES can leave
in its raw axis row. As for the mixed grid, `alpha=nfp*zeta`; physical toroidal
derivatives multiply the toroidal mode by `nfp` once.

There is no unique poloidal angle at `s=0`, so angular root solves still start
at `first_surface=1`. No fictitious Boozer axis ring or arbitrary poloidal
phase is serialized. cuMES supplies its native field on half surfaces rather
than at the axis. To let consumers compute their own field limit, v5 retains
the **unmodified first `radial_order` half surfaces**:

- `axis_half_s[j] = (j+0.5)/(source_ns-1)`, `j=0..radial_order-1`;
- `axis_half_sqrtg`, `axis_half_bsups`, `axis_half_bsupu`,
  `axis_half_bsupv`, `axis_half_bsubs`, `axis_half_bsubu`, `axis_half_bsubv`.

Each field has shape `[radial_order][source_nzeta][source_ntheta]`, with
native theta contiguous, at `theta=2*pi*k/source_ntheta`,
`alpha=2*pi*l/source_nzeta`. These are **native** fields, not remapped Boozer
fields. `u` and `v` denote native theta and physical zeta. Contravariant
components have units T/m, covariant components T m, and `sqrtg` has units
m^3 (the coordinates are dimensionless). In particular,

```text
|B|^2 = bsups*bsubs + bsupu*bsubu + bsupv*bsubv.
```

This stencil contains the data used for the transform's polynomial axis
extrapolation. A consumer can instead average over native theta before
extrapolating, or apply its own regularity treatment. Axis field values are
limits inferred from a finite radial grid; they are not exact extra samples.
The axis curve, signed profiles, native stencil, and non-axis Boozer fields
are all retained for that post-processing.

## Binary layout

Integers are little-endian signed int32; floating data are IEEE-754 binary64.
Strings are an int32 byte count followed by UTF-8 bytes. Arrays have no
length prefix except the explicit `aphi_count` below.

1. Eight bytes `MCBOOZ05`, then int32 version `5`.
2. Three strings: coordinate convention, Fourier convention, source path.
   The coordinate convention is identical to v3/v4. The Fourier convention
   remains `real-parity-v3` for symmetric and `real-parity-v4` for asymmetric.
3. Fourteen int32 values: `source_format_version`, `source_ns`,
   `source_ntheta`, `source_nzeta`, `source_mpol`, `source_ntor`, `nfp`,
   `first_surface`, `ntheta`, `nzeta`, `mmax`, `nmax`, `radial_order`, `lasym`.
   `lasym` must be 0 or 1.
4. Binary64 `resonance_tolerance`.
5. `s`, `iota`, each of length `surfaces=source_ns-first_surface`.
6. Int32 `mode_m`, `mode_n`, each of length `modes=(mmax+1)*(nmax+1)`.
7. Six arrays of length `surfaces*modes`: `rmncc`, `rmnss`, `zmnsc`,
   `zmncs`, `numnsc`, `numncs`. If `lasym=1`, append `rmnsc`, `rmncs`,
   `zmncc`, `zmnss`, `numncc`, `numnss` with the same extent.
8. `B`, `sqrt_g_b`, each of length `surfaces*nzeta*ntheta`, then
   `b2j00[surfaces]`.
9. String `tracing_convention` (the exact `TRACING_CONVENTION` in
   `src/boozer_output_internal.hpp`), int32 `aphi_count`, binary64 `phiedge`.
10. `aphi[aphi_count]`, `toroidal_flux_derivative[source_ns]`,
    `poloidal_flux_derivative[source_ns]`.
11. `axis_r_cos`, `axis_r_sin`, `axis_z_cos`, `axis_z_sin`, each of length
    `source_ntor+1`.
12. `axis_half_s[radial_order]`, then the seven native half-field arrays
    in the order listed above, each with
    `radial_order*source_nzeta*source_ntheta` entries.

## NetCDF and HDF5

Both use schema `magnetic-coordinate-boozer-v5`, integer attribute `lasym=0|1`,
and the existing coordinate/Fourier attributes. `phiedge` and
`tracing_convention` are additional root/global attributes. Arrays have the
names and extents above. Added NetCDF dimensions are `source_surface`,
`aphi_coefficient`, `axis_mode`, `axis_half_surface`, `source_zeta`, and
`source_theta`; HDF5 datasets have the equivalent shapes. All data remain real.

## Validation

Manufactured tests cover both field polarities and Jacobian orientations,
both radial orders, nonconstant flux profiles, asymmetric and symmetric axes,
and axis harmonics beyond the requested output truncation. They also verify
untouched native samples, polynomial axis limits, and binary/NetCDF/HDF5
serialization. Invalid dimensions, non-finite tracing data, inconsistent iota,
zero flux derivatives, noncanonical half grids, and truncated records are
rejected. Legacy v3/v4 round trips retain their missing-data state.
