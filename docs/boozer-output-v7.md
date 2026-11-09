# Boozer output format version 7

Version 7 retains all [v6 data](boozer-output-v6.md) and adds pressure and
the two covariant magnetic-field coefficients. The existing `iota` profile
remains in its original location. New transforms always provide these data.

## Surface profiles

Each array has shape `[source_ns-first_surface]` and shares the existing
`s` and `iota` grid. Its radial coordinate is `s=j/(source_ns-1)`; it is
proportional to toroidal flux only when `aphi` is constant.

| Name | Meaning | SI unit |
| --- | --- | --- |
| `pressure` | Prescribed equilibrium pressure | Pa |
| `current_i` | Poloidal covariant coefficient `I(s)` | T m |
| `current_g` | Toroidal covariant coefficient `G(s)` | T m |
| `iota` | Converged rotational transform, already present in v3-v6 | 1 |

The coefficients use the convention

```text
B = beta*grad(s) + I*grad(theta_b) + G*grad(zeta_b).
```

Here `zeta_b` is the physical toroidal angle. `I` and `G` are current
functions in T m; these arrays do not contain currents in amperes. No extra
`nfp`, `mu0`, or `2*pi` factor multiplies them.

The shift solver captures `I=<B_theta_native>` and
`G=<B_zeta_native>` from its normalized FFT zero modes, after radial
interpolation and any requested toroidal resampling. These are uniform
native-angle means, not averages on the mapped mixed grid. They retain the
magnetic-field sign and remain valid at rational iota. Field reversal
changes both signs and leaves pressure and iota unchanged.

The angular grid still omits the axis. Its curve, native field stencil,
and signed flux derivatives remain in the v6 tracing block.

## Pressure input

`CumesEquilibrium::pressure` and `CumesEquilibriumView::pressure` contain
`ns` finite values in Pa on the source full grid, including the axis.
Direct callers must provide this array; vacuum callers supply zeros.
The transform copies its non-axis entries without angular resampling.
Existing direct callers must populate this new input before transforming.

For a cuMES binary, the reader evaluates the embedded pressure prescription
on the full grid. With `gamma=adiabatic_index=0`, it uses

```text
t = min(s, spres_ped)
q = t * sum_k(aphi[k] * t^k)
x = min(abs(min(q, 1) * bloat), 1)
pressure(s) = pres_scale * P(x)
```

For `pmass_type=power_series`, `P(x)=sum_k(am[k]*x^k)`.
For `pmass_type=two_power`, `P(x)=am[0]*(1-x^am[1])^am[2]`, with at least
three coefficients. This matches cuMES's pressure prescription evaluated
at full-grid `s`; it does not interpolate its half-grid pressure samples.
The profile argument is not divided by `q(1)`. Pressure uses Pa, without
the solver's internal `mu0` factor.

The reader rejects unsupported profile types, nonzero gamma, incomplete
`two_power` coefficients, and non-finite inputs or intermediate values.
The transform checks the pressure dimension and finiteness. Output
validation also checks the dimensions and finiteness of all three new
arrays and rejects zero or cancelled `G+iota*I`.

## Containers and compatibility

Binary output uses magic `MCBOOZ07` and int32 version `7`. The rest of the
v6 layout is unchanged. After its final `nu_grid` array, v7 appends:

1. The length-prefixed UTF-8 profile convention below.
2. `pressure`, as `source_ns-first_surface` binary64 values.
3. `current_i`, with the same length and type.
4. `current_g`, with the same length and type.

```text
surface-profiles-v1: pressure in Pa; B=beta*grad(s)+I*grad(theta_b)+G*grad(zeta_b); current_i=I and current_g=G in T m with physical zeta_b; all profiles share the non-axis s/iota grid
```

NetCDF and HDF5 use schema `magnetic-coordinate-boozer-v7`, a
`profile_convention` attribute, and real one-dimensional arrays with the
names above. Their `units` attributes are `Pa`, `T m`, and `T m`.
`iota` keeps its existing dimensionless representation.

In memory, `BoozerResult::profiles` holds `BoozerSurfaceProfiles`.
Readers accept v7/v8 and can select pressure and currents independently.
Writers emit v8 and require profiles, tracing, magnetic fields, and maps.

Tests cover analytic pressure profiles, both pressure types, pedestal and
flux clipping, both radial orders, field reversal, nonunit `nfp`, rational
iota, HOST/CUDA shift coefficients, all three containers, optional quantity failures,
and malformed or truncated profile data.
