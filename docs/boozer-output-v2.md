# Boozer output format version 2

The binary, NetCDF, and HDF5 containers implement the same version-2 data
model. Files contain only signed integers, UTF-8 metadata, and real IEEE-754
binary64 values. Complex numbers are not part of the public format; the
complex coefficients used by the internal FFT are folded into real Fourier
parity families before writing.

The coordinate convention is the mixed grid

```text
(s, theta_b, zeta),       zeta_b = zeta + nu(theta_b, zeta),
```

where `theta_b` is uniform and `zeta` is the unchanged source toroidal angle.
The magnetic axis is excluded (`first_surface=1`) because flux angles are
degenerate there.

## Real Fourier representation

Modes are nonnegative and n-major: `n=0..nmax`, then `m=0..mmax`. The
toroidal integer is a field-period mode. The stored real amplitudes reconstruct
the mixed-grid fields as

```text
R  = sum [rmncc cos(m theta_b) cos(n zeta)
        + rmnss sin(m theta_b) sin(n zeta)]

Z  = sum [zmnsc sin(m theta_b) cos(n zeta)
        + zmncs cos(m theta_b) sin(n zeta)]

nu = sum [numnsc sin(m theta_b) cos(n zeta)
        + numncs cos(m theta_b) sin(n zeta)].
```

The transform is stellarator-symmetric: `R` is even, while `Z` and `nu` are
odd. A result whose internal spectrum does not satisfy these real parity
relations is rejected instead of serializing a lossy representation. FFTs are
unnormalized internally, but the six stored families are reconstruction
amplitudes; consumers do not apply an FFT normalization factor.

## Common data model

Every container records:

- coordinate, Fourier, source, grid, truncation, interpolation, and resonance
  metadata;
- `s[surface]`, `iota[surface]`, and `b2j00[surface]`;
- `mode_m[mode]`, `mode_n[mode]`, and the six real coefficient arrays
  `[surface][mode]`;
- real `B[surface][zeta][theta_b]` and
  `sqrt_g_b[surface][zeta][theta_b]` arrays.

The last index is contiguous. The zero mode and Jacobian obey

```text
b2j00 = integral B^2 sqrt(g_p) dtheta_p dzeta,
sqrt_g_b = b2j00 / (4*pi^2*B^2).
```

Thus `b2j00` is the unnormalized integral, not the flux average. The sign of
`sqrt_g_b` preserves the source Jacobian orientation.

## Binary container

The little-endian binary begins with `MCBOOZ02` and a signed 32-bit version
equal to 2. Strings are a signed 32-bit byte count followed by UTF-8 bytes.
After the metadata header, the payload order is `s`, `iota`, `mode_m`,
`mode_n`, the six coefficient families in the order shown above, `B`,
`sqrt_g_b`, and `b2j00`.

## NetCDF and HDF5 containers

NetCDF dimensions are `surface`, `theta_b`, `zeta`, and `mode`. HDF5 datasets
have the equivalent extents. Both use the common names above and root/global
attributes for metadata. Coefficients are ordinary rank-2 double arrays; no
compound datatype or real/imaginary component dimension is used.
