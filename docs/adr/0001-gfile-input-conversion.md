# ADR-0001: G-file conversion to cuMES input

Date: 2026-11-18

## Decision

Add a separate optional `gfile` component. Keep libmeq's `delegation` branch
as a pinned submodule and consume its headers privately. Its CMake project
changes compiler flags and builds tests unconditionally, so the adapter does
not add that project to the build. Public types use the existing normalized
`InputParams` and plain LCFS points; they expose no libmeq or spline types.
This keeps the dependency direction from cuMES to equilibrium-io.

Use libmeq's corrected COCOS transformations, reader, contour finder, and
spline dependency. A bounded preflight guards assumptions in its reader.
Map signs and flux units to COCOS 5, matching cuMES's angles and negative
Jacobian. Retain both physical field and current polarities.

Integrate the piecewise linear q profile and remap pressure and iota to uniform
normalized toroidal flux. Use fixed-iota input with `aphi={1}`. Fit Fourier
boundary modes after periodic radius interpolation. Preserve the supplied
LCFS by default; an explicit flux crop must bracket every contour ray and
pass a residual check. Report sampled boundary-fit and integrated
poloidal-flux errors, and reject results outside the requested tolerances.

Measure LCFS shape separately on the original polygon. Use its radial extrema
for center and minor radius, and magnetic-axis height as the default upper/lower
reference. Do not impose conversion-only profile or star-shapedness constraints
on geometry reading and measurement.

## Qualification

This is a Class C input-conversion algorithm. It leaves solver arithmetic and
existing I/O formats unchanged. Qualification uses independent analytic files,
all 64 COCOS/polarity combinations, interpolation and Fourier refinement,
flux identities, malformed-input checks, sanitizer tests, and cuMES's actual
parser and radial evaluator. These checks establish the conversion contract;
they do not establish equilibrium convergence or browser execution.

See [the G-file guide](../gfile.md) for formulas, units, bounds, and API examples.
