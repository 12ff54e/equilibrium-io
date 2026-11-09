# Boozer output format version 6

> Historical schema: current readers reject this file version. This document
> defines scientific blocks retained by v7/v8. Current API and admission rules
> are in [File compatibility](compatibility.md).

Version 6 exports all six magnetic-field components in the **Boozer basis**
`(s, theta_b, zeta_b)`, together with their signed Jacobian and the untruncated
coordinate map. Magnetic-coordinate performs the basis transformation;
consumers such as GEPT perform interpolation, derivatives, and particle
tracing. All [v5 data](boozer-output-v5.md), including signed flux profiles,
the converged axis, and its native field stencil, remain present.

[Version 7](boozer-output-v7.md) adds pressure and the covariant coefficients
`I` and `G` to this layout. New transforms write v7. Binary readers continue
to accept v3/v4/v5/v6, with missing data explicitly absent.

## Components and sampling

The following arrays have shape `[source_ns-first_surface][nzeta][ntheta]`,
with Boozer theta contiguous. They share the non-axis mixed grid used for
`b`, `sqrt_g_b`, and the geometry. In particular, **the component basis and
the grid's independent toroidal variable are different**: components use
physical `zeta_b`, while samples retain `alpha=nfp*zeta_native`.

| Name | Meaning | SI unit |
| --- | --- | --- |
| `boozer_sqrtg` | Signed Jacobian for `(s,theta_b,zeta_b)` obtained from the coordinate map | m³ |
| `boozer_bsups` | `B dot grad(s)` | T/m |
| `boozer_bsupu` | `B dot grad(theta_b)` | T/m |
| `boozer_bsupv` | `B dot grad(zeta_b)` | T/m |
| `boozer_bsubs` | `B dot partial(x)/partial(s)` at fixed Boozer angles | T m |
| `boozer_bsubu` | `B dot partial(x)/partial(theta_b)` at fixed `(s,zeta_b)` | T m |
| `boozer_bsubv` | `B dot partial(x)/partial(zeta_b)` at fixed `(s,theta_b)` | T m |

Angles and source radial `s=j/(source_ns-1)` are dimensionless. Do not replace
`s` with toroidal flux without the corresponding radial change of basis.
In particular, `B^s=0` does not imply `B_s=0`; the radial covariant component
is retained. The ideal-MHD conversion assumes zero radial current
`j dot grad(s)`, which also does not require `B_s=0`.
Components change sign under field reversal, while the signed
Jacobian does not. Neither a positive orientation nor stellarator symmetry
is imposed on this block.

The untruncated maps have the same shape as the components:

- `source_theta_shift = theta_native - theta_b` is an unwrapped displacement.
- `nu_grid = zeta_b - zeta_native` is in physical toroidal radians.

At each sample,

```text
theta_native = theta_b + source_theta_shift   (mod 2*pi),
zeta_native = alpha/nfp,
zeta_b = alpha/nfp + nu_grid.
```

The poloidal displacement is computed as `-lambda(theta_native)-iota*nu`,
without subtracting wrapped angles. This preserves smooth derivatives even
when either angle crosses its periodic boundary. The existing truncated
Fourier representation of `nu` remains available separately.

## Basis transformation

Native half-grid components and the signed native Jacobian are first
unstaggered with the selected two- or four-point radial interpolation.
Existing periodic cubic B-spline interpolation resamples source zeta and
then evaluates native theta at the mapped output positions.

Let `p=(s,u,v)=(s,theta_b,zeta_native)`, `d=source_theta_shift`, and
`nu=nu_grid`. Differentiate these smooth maps on the mixed grid. Map
interpolation is cubic where at least four samples are available, quadratic
for three samples, and linear for two; a singleton toroidal dimension has
zero derivative. Radial derivatives use the non-axis surfaces only, with
one-sided endpoint information. Subscripts below mean partial derivatives
holding the other mixed coordinates fixed. Physical `v` derivatives apply
`nfp` once to derivatives with respect to `alpha`.

For native components `Q`, first transform from `q=(s,u+d,v)` to `p`:

```text
M^s = Q^s,
M^u = (Q^theta - d_s*Q^s - d_v*Q^zeta)/(1+d_u),
M^v = Q^zeta,

M_s = Q_s + d_s*Q_theta,
M_u = (1+d_u)*Q_theta,
M_v = Q_zeta + d_v*Q_theta.
```

Then transform to `b=(s,u,v+nu)`:

```text
B^s      = M^s,
B^theta_b = M^u,
B^zeta_b  = nu_s*M^s + nu_u*M^u + (1+nu_v)*M^v,

B_s       = M_s - nu_s*M_v/(1+nu_v),
B_theta_b = M_u - nu_u*M_v/(1+nu_v),
B_zeta_b  = M_v/(1+nu_v),

boozer_sqrtg = J_native*(1+d_u)/(1+nu_v).
```

Non-finite maps, non-positive map determinants at output samples, vanishing
Jacobians, and invalid transformed field strengths are rejected. The basis
transformation preserves `B^i B_i` of the sampled native components. It does
not average angular covariant components into ideal flux functions or
renormalize the field to enforce an identity.

At finite resolution, separately interpolating the components, `B^2`, and
Jacobian is not algebraically identical to interpolating their products.
Consequently `sum(boozer_B^i*boozer_B_i)` can differ from the existing `b^2`,
and `boozer_sqrtg` can differ from legacy `sqrt_g_b=b2j00/(4*pi^2*b^2)`.
The latter retains its original normalization and byte values. For a
consistent tensor calculation use the component block with `boozer_sqrtg`;
assess interpolation and angular/radial resolution for the intended orbit.
The output does not assert exact force balance or tracing accuracy.

The angle conventions agree with the
[Boozer transformation theory](https://hiddensymmetries.github.io/booz_xform/theory.html).

## Guiding-center derivatives

The saved `boozer_bsupu` directly provides `B dot grad(theta_b)`. For any
scalar `f` sampled on the mixed grid, including `|B|`, its derivatives at
fixed Boozer coordinates are

```text
D_s f       = f_s - nu_s*f_v/(1+nu_v),
D_theta_b f = f_u - nu_u*f_v/(1+nu_v),
D_zeta_b f  = f_v/(1+nu_v).
```

Using coordinate order `(s,theta_b,zeta_b)` and the signed Jacobian,

```text
(B cross grad(|B|)) dot grad(theta_b)
    = (B_zeta_b*D_s|B| - B_s*D_zeta_b|B|)/boozer_sqrtg.
```

Consumers can instead remap every field onto a uniform true Boozer toroidal
grid before differentiating. Derivatives at fixed source `alpha` must not
be treated as derivatives at fixed `zeta_b`.

The axis remains a separate curve, outside the angular root solves and
component transformation. The v5 `axis_half_*` arrays still contain native
near-axis samples for the consumer's axis limit; they are explicitly **not**
Boozer-basis components on an axis ring.

## Independent consistency checks

Reconstruct the field strength from all six components:

```text
B2_components = B^s*B_s + B^theta_b*B_theta_b + B^zeta_b*B_zeta_b.
```

The original toroidal grid is also sufficient to reconstruct the signed
Boozer Jacobian from the exported `R,Z` geometry. With mixed coordinates
`(s,u,v)=(s,theta_b,zeta_native)`, the cylindrical position is
`x=(R*cos(v),R*sin(v),Z)`. Its volume determinant is

```text
J_mixed = R*(R_u*Z_s - R_s*Z_u).
```

Since `zeta_b=v+nu(s,u,v)`, the determinant of the remaining coordinate
change is `1+nu_v`. Therefore

```text
J_geometry = R*(R_u*Z_s - R_s*Z_u)/(1+nu_v),
nu_v = nfp*partial(nu_grid)/partial(alpha).
```

All derivatives in this expression hold the other mixed coordinates fixed.
Neither `R_v` nor `Z_v` is needed. No remapping to a uniform Boozer toroidal
grid is required. In axisymmetry `nu_v=0`, so the uncorrected cylindrical
determinant already equals the Boozer Jacobian; in general 3-D it does not.

The identity `B^2 J=F(s)` is required of every valid Boozer map, including
generalized Boozer coordinates with nonzero radial current. The
[covariant construction](https://github.com/12ff54e/magnetic-coordinate/blob/main/docs/adr/0001-covariant-boozer-shift.md) assumes zero
radial current, but this check applies independently of that assumption.

Compare `B2_components` with `b^2`, and `J_geometry` with both
`boozer_sqrtg` and the separately normalized `sqrt_g_b`. The independent
product `B2_components*J_geometry` should approach
`b2j00[surface]/(4*pi^2)` with sufficient radial/angular resolution and
retained Fourier harmonics. Checking only stored `b^2*sqrt_g_b` cannot expose
a conversion defect, because that identity is imposed by construction.
Truncating the transformed geometry at the source mode cutoff adds error:
the nonlinear angle map can generate additional harmonics needed for an
accurate geometry derivative.

The production `transform_to_boozer` path calls
`check_boozer_consistency` before returning a result on HOST, CUDA, and Wasm.
It rejects non-finite or zero geometry Jacobians, orientation mismatches with
the component Jacobian, and folded toroidal maps. The geometry uses the
retained output spectrum, so truncation cannot silently discard information
needed to preserve orientation. The error identifies the source surface and
angular sample; a failed browser conversion leaves VMEC available and Boozer
disabled. Existing binary files can be checked with the same public function.

The checker reports maximum and RMS relative errors for `B^2`, both
Jacobians, and `B^2 J`. It does not impose the manufactured-case accuracy
bounds on unrelated equilibria. Passing orientation alone does not establish
a valid Boozer map or orbit accuracy; a `B^2 J` discrepancy remains a violation
of the defining coordinate condition. No field or geometry is renormalized to
pass the check.

The [near-rational conditioning diagnosis](https://github.com/12ff54e/magnetic-coordinate/blob/main/docs/boozer-conditioning.md) explains
the former W7-X shift failure and distinguishes it from FFT backend
differences and retained-geometry truncation. The
[covariant-field construction](https://github.com/12ff54e/magnetic-coordinate/blob/main/docs/adr/0001-covariant-boozer-shift.md) now avoids
that rational-surface divisor. File layout, coordinate basis, zero-mode
normalization and PEST mean-zero gauge are unchanged. The
[direct native-to-Boozer remap](https://github.com/12ff54e/magnetic-coordinate/blob/main/docs/adr/0002-direct-native-boozer-map.md) evaluates
that gauge as `mean_native(nu*(1+lambda_theta))=0` and the zero mode as
`4*pi^2*mean_native(B^2 J_native)`, without a uniform PEST grid. These preserve
the continuous definitions; discrete values differ from the former two-stage
spline resampling. The legacy
`resonance_tolerance` field is retained as metadata and is unused by the new
construction; readers cannot infer an algorithm revision from the file version.

`magnetic-coordinate.geometry-consistency` runs the entire transform and
binary round trip on analytic toroidal geometry with metric-consistent,
divergence-free input fields. These manufactured fields test tensor/geometry
transformation. Coverage includes axisymmetry, symmetric and asymmetric 3-D,
signed orientations and field polarities, and refinement
from 17 to 33 full radial surfaces and 32 to 64 poloidal points. Three-D
cases use `nfp=3` and 16 to 32 toroidal points. Geometry is synthesized from
the serialized coefficients; Fourier collocation derivatives and a local
six-point polynomial in `sqrt(s)` provide independent angular and radial
checks, including the first non-axis surface and LCFS.

The geometry regression requires pointwise relative errors below `6e-5`
for the component contraction and `1e-3` for the geometry/component Jacobian,
plus decreasing maximum and RMS errors under refinement. Its general
divergence-free fields have nonzero radial current and a nonclosed angular
covariant one-form. The reported `B^2 J` errors identify this construction's
failure to produce Boozer maps for those inputs; the errors are not permitted
departures from the Boozer condition. Compatible analytic pullbacks of
prescribed smooth Boozer maps provide positive tests of that identity,
constant angular covariant components, and straight field-line pitch in
`magnetic-coordinate.host` and `magnetic-coordinate.boozer-fft`, including
exact rational iota. The geometry regression also detects deliberately
corrupted components, geometry, and a missing toroidal map. These bounds
qualify the specified manufactured cases, not arbitrary equilibria or orbit
accuracy.

The same executable also converts and serializes an exact vacuum torus
`B=C/R e_zeta` with an analytic nonlinear PEST map and zero iota. Its
symmetric and poloidally rotated asymmetric representations check contraction
to `2e-13`, both geometry Jacobians and `B^2 J` to `1e-3`, the known flux
constant, and improvement under radial/angular refinement. This retains an
end-to-end geometry/field Boozer-identity check on a physical equilibrium
satisfying the covariant construction's assumptions.

With testing enabled, the same executable reports relative maximum/RMS
errors for existing binary v6 files:

```bash
./build/magnetic_coordinate_geometry_consistency_test boozer-output.bin
```

External-file mode reports diagnostics without applying manufactured-case
acceptance thresholds. It rejects invalid data, folded toroidal maps, and
opposite geometry/component Jacobian orientations.

## Serialization and compatibility

The public API adds `BoozerTracingData::magnetic`, containing
`BoozerMagneticData::fields` in the seven-field table order. Sampled maps are
`BoozerMixedGrid::source_theta_shift` and `BoozerMixedGrid::nu`. The optional
magnetic block is absent in legacy v3/v4/v5 binary files.

Binary v6 starts with `MCBOOZ06` and int32 version `6`. The header and payload
through the v5 axis stencil are unchanged apart from the magic/version.
Append, in this order:

1. String `magnetic_convention`, equal to `MAGNETIC_CONVENTION` in
   `src/boozer_output_internal.hpp`.
2. Seven arrays in the table order, each containing
   `(source_ns-first_surface)*nzeta*ntheta` doubles.
3. `source_theta_shift` and `nu_grid`, each with the same extent.

Strings use the existing int32 byte-length plus UTF-8 bytes. All new numeric
values are little-endian IEEE-754 binary64. There is no native half-grid
component block beyond the retained v5 axis stencil.

NetCDF and HDF5 use schema `magnetic-coordinate-boozer-v6` and root/global
attribute `magnetic_convention`. The new variables/datasets use the table
names and existing dimensions `surface,alpha,theta_b`.

Rewriting a result with tracing but without `magnetic` writes v5; a result
without tracing retains v3/v4. Readers never invent missing components.
Tests cover both symmetries and all three containers, malformed/truncated
payloads, field reversal, both Jacobian orientations, nonzero radial
covariant fields, unwrapped displacements, endpoint map derivatives, and
independent Cartesian checks of the components and guiding-center drift
contraction.
