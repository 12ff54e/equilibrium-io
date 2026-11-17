# G-file input conversion and LCFS shape

The optional `equilibrium::gfile` component reads GEQDSK files with libmeq and
creates fixed-boundary, fixed-iota cuMES input. It also measures the upper and
lower LCFS shape. It does not solve an equilibrium or certify convergence.
This feature uses the pre-0.4 equilibrium-io API on the `miller-batch` branch.

## Build and use

The `deps/libmeq` submodule tracks `delegation` and is pinned to a commit that
includes its COCOS sign and flux-unit fixes. Initialize its spline dependency:

```sh
git submodule update --init --recursive deps/libmeq
cmake -S . -B build -DEQUILIBRIUM_USE_LIBMEQ=ON
cmake --build build
```

`EQUILIBRIUM_USE_LIBMEQ` accepts `AUTO` (default), `ON`, or `OFF`. `AUTO`
builds the component when both submodules are present; `ON` requires them.
The implementation consumes libmeq headers privately. Installed users need
only the compiled library, its public headers, and a C++20 standard library.
Selecting `data`, `binary`, or `io` does not load `gfile`. No optional codecs
or CUDA libraries are needed by `gfile`. Installations include both dependency
license notices under `share/equilibrium_io/licenses`.

```cmake
find_package(equilibrium_io CONFIG REQUIRED COMPONENTS gfile)
target_link_libraries(my_converter PRIVATE equilibrium::gfile)
```

The standalone tool writes cuMES configuration JSON or prints shape metrics:

```sh
equilibrium-gfile g123456 --cocos 1 --output input.json
equilibrium-gfile g123456 --cocos 1 --psi-ratio 0.99 --output inner.json
equilibrium-gfile g123456 --shape --cocos 1
```

`--shape` writes JSON to standard output and cannot be combined with `--output`.
Other conversion options are `--mpol`, `--ntheta`, `--profile-points`, and
`--boundary-points`. Use the C++ settings to change fit tolerances or the solver
stage. Default settings are `mpol=12`, `ntheta=30`, 512 boundary samples, 129
profile samples, and one stage with 65 surfaces, 10000 iterations, and `1e-10`
tolerance. Custom stages require 3–512 surfaces and a tolerance of at least
`1e-16`; resolved `ntheta` must be even and at most 256. These are starting settings for a native double solve, not evidence
that a given equilibrium converges.

```cpp
#include <equilibrium/gfile.hpp>

equilibrium::GfileInputSettings settings;
settings.cocos = 1;
settings.psi_ratio = 0.99;
auto converted = equilibrium::gfile_to_cumes_input("g123456", settings);
const equilibrium::InputParams& input = converted.input;
std::string json = equilibrium::gfile_input_to_json(converted);
auto status = equilibrium::write_gfile_input("input.json", converted);
if (!status) throw std::runtime_error(status.error());
auto shape = equilibrium::calculate_lcfs_shape(converted.lcfs);
```

File and stream overloads share the same conversion. Reads and calculations
throw on failure; the atomic writer returns `Status`. The returned `InputParams`
has resolved dimensions, axis arrays, and matching raw/folded Fourier modes.
The JSON uses the flat cuMES input schema, rather than a scientific-result
container. Existing output is preserved if validation or writing fails before
publication.

## Coordinate and profile contract

libmeq determines source COCOS from an explicit API override, a header tag,
a producer hint, or sign inference, in that order. Inference is heuristic:
unknown producers default to EFIT-style sign rules and cannot identify an
unlabelled extra factor of 2π. Supply the known COCOS index when possible.
Indices 1–8 store poloidal flux per radian; 11–18 store full flux.

The converter uses COCOS 5: toroidal angle is counterclockwise from above,
and poloidal angle is counterclockwise in the R–Z plane. With signed
`psi_5` in Wb/radian, `Delta_psi_5=psi_edge-psi_axis`, and source normalized
poloidal flux `u`, it uses

```text
phiedge = -2*pi*Delta_psi_5 * integral_0^psi_ratio q_5(u) du
iota    = 1/q_5
```

This follows cuMES's negative coordinate Jacobian and
`Phi'(s)=-phiedge/(2*pi)`, `chi'(s)=iota(s)*Phi'(s)` when `aphi={1}`.
Both independent current and toroidal-field polarities are retained.
`sign(phiedge)` must agree with the converted `F=R*B_phi`. The sign definitions
follow [Sauter and Medvedev, Table I](https://crppwww.epfl.ch/~sauter/cocos/Sauter_COCOS_Tokamak_Coordinate_Conventions.pdf).

The input uses `ncurr=0`. Its pressure and iota lists are sampled uniformly in
normalized **toroidal** flux, not at the G-file's poloidal-flux samples.
The converter integrates the piecewise linear q profile exactly and inverts
that monotone map, then linearly interpolates pressure and q before taking
`1/q`. Pressure stays in Pa, `pres_scale=1`, and the stored current is converted
to the same toroidal orientation but is unused by fixed-iota mode.

The output list introduces interpolation error. The converter independently
integrates the output iota list and checks the resulting edge poloidal flux
against `psi_ratio*Delta_psi_5`. `poloidal_flux_relative_error` reports this
relative discrepancy; the default limit is `1e-3`. Increase `profile_points`
when needed. No q zero, q-sign reversal, inconsistent F sign, or negative
pressure is accepted for input conversion.

## Boundary fitting and cropping

At `psi_ratio=1`, the supplied LCFS is used without a flux-grid contour search.
For a smaller ratio, libmeq finds intersections along rays from the magnetic
axis through LCFS vertices. Each ray must bracket the requested flux, and each
returned point must pass a flux residual check. Grid-domain bounds are checked;
the requested crop is never silently changed to fit an inconsistent file.
Cropping requires a grid of at least 6×6 points.

The boundary must be ordered and star-shaped about the magnetic axis. A
periodic cubic spline of geometric radius provides uniform angular samples;
Fourier modes `0 <= m < mpol` then describe R and Z. Absolute axis height and
up/down asymmetry are retained through `lasym`, `rbs`, `zbc`, and `zaxis_c`.
cuMES supports `2 <= mpol <= 16`.

`boundary_max_error` is the largest distance in metres between the Fourier
curve and the sampled contour. Its default limit is `1e-3` times the original
LCFS minor radius. Positive R and boundary orientation are checked on the
sampling grid. These are sampled checks, not error bounds between samples;
source resolution, interpolation, and Fourier truncation each affect accuracy.
A sharp separatrix may need an inner crop to fit within the mode limit.

## Upper and lower LCFS shape

Shape uses the original supplied polygon, even when conversion crops the
boundary. It requires neither a valid q/F profile nor a star-shaped contour.
Repeated closing vertices and consecutive identical points are omitted.
No smoothing, cropping, or Fourier truncation enters the measurements.

```cpp
auto lcfs = equilibrium::read_gfile_lcfs("g123456");
auto shape = equilibrium::calculate_lcfs_shape(lcfs);
// Or choose a different horizontal reference explicitly:
auto custom = equilibrium::calculate_lcfs_shape(lcfs.points, z_reference);
```

The default horizontal reference is the **magnetic-axis height** `Z_axis`.
The radial center is geometric, from the LCFS's global radial extrema:

```text
R_center = (R_max + R_min)/2       a = (R_max - R_min)/2
kappa_upper = (Z_max - Z_axis)/a  kappa_lower = (Z_axis - Z_min)/a
delta_upper = (R_center - R_at_Z_max)/a
delta_lower = (R_center - R_at_Z_min)/a
```

Positive triangularity means that the upper or lower tip lies inward of
`R_center`; negative triangularity is retained. If several vertices tie at
an extreme height, their minimum and maximum R define a midpoint. This makes
flat-top and flat-bottom results independent of point density. The chosen
Z reference must lie strictly between the vertical extrema. Dimensions are in
metres; elongation and triangularity are dimensionless. This definition does
not use the magnetic-axis R, equatorial crossings, or an average of both halves.

## Validation and limits

Tests cover all 16 COCOS conventions with four field/current polarity pairs,
profile and Fourier refinement, crop failures, the output flux integral,
shape reflection/translation, unequal half-shapes, flat extrema, repeated
vertices, reversed order, and malformed input. Installed consumers exercise
the public API without libmeq headers. Generated symmetric and asymmetric
JSON is also checked with cuMES's actual strict parser and radial evaluator.

The reader bounds input to 128 MiB, grid dimensions to 2–4096 with at most
2^22 grid points, and boundary/limiter counts to 100000. It checks complete finite
numeric data before invoking libmeq's reader and accepts Fortran `D` exponents.
It expects the standard 48-character GEQDSK header followed by three I4 fields.
Additional producer data after the standard payload are not interpreted.
