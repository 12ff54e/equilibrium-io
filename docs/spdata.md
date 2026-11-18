# Axisymmetric spdata output

`equilibrium/spdata.hpp` defines the owning `Spdata` record and borrowed
`SpdataView`. `equilibrium/spdata_validation.hpp` provides header-only
validation. Link `equilibrium::io` for `write_spdata(path, record)` from
`equilibrium/spdata_output.hpp`. No optional codec is required.

The writer serializes prepared quadratic coefficients. The producer owns
coordinate conversion, profile calculation, spline fitting, units, and
physical validity. `magnetic-coordinate` is the intended producer. This is
an axisymmetric format; the record has no toroidal dimension or ripple data.
It imposes no up/down symmetry.

```cpp
#include <equilibrium/spdata_output.hpp>

// Fill all coefficients in the units required by the receiving application.
equilibrium::Spdata data = prepared_coefficients;
auto status = equilibrium::write_spdata("spdata.dat", data);
if (!status) throw std::runtime_error(status.error());
```

`as_view(data)` borrows its arrays and metadata; the owner must remain alive
and unchanged during use. Calls with temporary owners are disabled. A caller
with existing arrays may fill `SpdataView` directly.

## Grid and coefficients

The layout follows MCT's `Spdata::print` and `SpdataLiteral` in `Spdata.h`
and `Spdata.cpp`:

- `lsp` is the number of radial records, including the axis. Its allowed
  range is 2 through 9999. `delta_psi = psi_wall / (lsp - 1)`.
- `lst` is the number of unique poloidal points, from 1 through 999.
  `delta_theta = 2*pi/lst`. The writer appends a copy of the first point to
  every coefficient row. It never stores that endpoint in the record.
- `psi_wall` and `psi_sep` are separate, finite, positive flux limits. The
  writer preserves both; it does not assume they are equal.
- `b`, `r`, `z`, and `j` hold field strength, cylindrical R, vertical Z, and
  the Jacobian in the producer's poloidal-flux/Boozer-angle coordinates.
  Each contains `lsp * 9 * lst` values in `[radial][coefficient][theta]`
  order. `field_index()` gives that flat index.
- `q`, `f`, `i`, `p`, `r_minor`, and `psi_t` hold safety factor, poloidal
  current function, toroidal current function, pressure, normalized outer
  minor radius, and toroidal flux. Each contains `lsp * 3` values in
  `[radial][coefficient]` order. `profile_index()` gives that flat index.

For each 2-D cell, coefficient `a + 3*b` multiplies
`x^a * delta_theta_local^b`, with `a,b` from zero through two. These are
polynomial coefficients, **not raw derivatives**: a second derivative
contributes half its value. Here `x = sqrt(psi)` in the first radial record;
later records use `x = psi - radial*delta_psi`. The first axis value is
constant in theta, so coefficients 3 and 6 must be zero there.

The 1-D coefficients multiply `1, x, x*x`. Only `r_minor` uses `sqrt(psi)`
in its first record. The remaining 1-D fields use `psi` there. All later
records use the same local `delta_psi` as the 2-D fields. Producers must
supply the final radial record too, even though it starts at `psi_wall`.
The writer neither refits nor rescales any coefficient.

## Text contract

The first line is `meta_data` (which may be empty). Embedded CR, LF, and NUL
are rejected, as are metadata strings longer than 1 MiB. The next line has
four I4 integers: `lsp`, `lst`, `4`, `8`. The following line holds `psi_wall`
and `psi_sep` as uppercase scientific E18.10 values.

Each radial record writes `b,r,z,j` in that order. For each field it writes
nine coefficient rows, each with `lst+1` values, wrapping after four values.
It then writes `q,f,i,p,r_minor,psi_t`, three coefficients per line. All real
values use E18.10 and a period as the decimal separator regardless of locale.
Decimal output rounds doubles to eleven significant digits.

The footer contains two I4 zeros, then `R_axis,0,0`, then `0,0`, using the
same real format. `R_axis` comes from the first constant R coefficient and
must be positive. The output contains no ripple. MCT's literal reader skips
the footer; this writer retains it for downstream consumers.

The `lst <= 999` limit preserves compatibility with MCT's whitespace parser:
a four-digit second integer would join the first integer in its I4 header.
The first integer may use all four columns. No delimiter extension is used.

## Units and publication

The record deliberately carries no implicit unit conversion. MCT's usual
normalized output uses axis R and B as scales; it shifts Z by the axis
height, normalizes the outer minor radius by axis R, and writes pressure
as `mu0*p/B0^2`. Its direct physical-unit path writes pressure as `mu0*p`,
whereas its literal `convert_to_SI` helper leaves pressure unchanged.
Producers must choose and describe one convention; the writer cannot infer
it from coefficients. `meta_data` can record that choice.

Validation checks dimensions, shapes, finite values, axis constants, metadata,
and positive flux limits and axis radius. It does not establish force balance,
polynomial continuity, or convergence. Validation runs before opening output.
The writer uses a same-directory temporary file and the shared checked
flush/fsync/rename protocol. Failures before rename preserve the destination
and remove the temporary file. A directory-sync failure after rename reports
that replacement occurred but durability could not be confirmed.

Use this dedicated API for spdata. Native output suffix dispatch and
`equilibrium-convert` continue to handle complete native/Boozer records;
they do not perform the numerical conversion needed to form these coefficients.
