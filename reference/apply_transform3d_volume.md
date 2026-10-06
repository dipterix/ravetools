# Apply a registration (linear and/or deformable) to a volume or to points

`apply_transform3d_volume` resamples a 3D volume (or a 4D stack of
frames, each frame treated identically) through a chain of transforms
onto a reference grid; `apply_transform3d_points` maps `RAS` coordinates
through the same chain. The chain is either a registration result (from
[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md)
or
[`load_registration`](https://dipterix.org/ravetools/reference/save_registration.md)),
whose mapping direction is chosen with `direction`, or an 'ANTs'-style
list of transforms applied in the given order. Both functions stream
every voxel or point through the chain in `C++`; no composite
deformation field is ever materialized. These functions supersede
[`apply_transform3d`](https://dipterix.org/ravetools/reference/apply_transform3d.md),
which only handles a single linear transform.

## Usage

``` r
apply_transform3d_volume(
  volume,
  transforms,
  vox2ras = NULL,
  reference_dim = NULL,
  reference_vox2ras = NULL,
  direction = c("forward", "inverse"),
  invert = FALSE,
  interpolation = c("trilinear", "nearest"),
  na_fill = 0
)

apply_transform3d_points(
  points,
  transforms,
  direction = c("forward", "inverse"),
  invert = FALSE
)
```

## Arguments

- volume:

  a 3D array to resample (integer and logical arrays are converted to
  double), or a 4D array whose frames are all transformed the same way;
  trailing dimensions of size one beyond the fourth are dropped

- transforms:

  a `ravetools_register_volume3d` object, or a list of transforms in
  'ANTs' order (see the sections below)

- vox2ras:

  the \\4\times 4\\ voxel-to-`RAS` matrix of `volume`; default `NULL`
  looks for the array's `"vox2ras"` attribute, then at the registration
  geometry. This and every other matrix involved (reference grid, field
  grids, transforms) must be finite and non-singular; a singular one is
  an error

- reference_dim, reference_vox2ras:

  the output grid: its dimension (length 3) and voxel-to-`RAS` matrix;
  required for a transform list, defaulting to the registration's target
  (`"forward"`) or source (`"inverse"`) grid otherwise

- direction:

  which mapping of a registration object to evaluate (see the table);
  only valid with a registration object

- invert:

  logical, recycled to the list length: invert the corresponding matrix
  before use; only valid with a transform list, and only for matrices

- interpolation:

  `'trilinear'` (default) or `'nearest'` (for labels or masks)

- na_fill:

  value for output voxels that fall outside `volume`; default `0`

- points:

  an `N x 3` matrix, a data frame whose first three columns are
  `x, y, z`, or a length-3 vector, in `RAS` millimeters; rows with a
  missing or non-finite coordinate give `NA` rows

## Value

`apply_transform3d_volume` returns the resampled volume with dimension
`reference_dim` (plus the frame dimension for 4D input) and a
`"vox2ras"` attribute equal to `reference_vox2ras`.
`apply_transform3d_points` returns an `N x 3` double matrix.

## Registration objects

A `ravetools_register_volume3d` object holds the linear transform \\A\\
(target `RAS` to source `RAS`) and, for `"syn"` registrations, the
forward and inverse displacement fields \\u_f\\ and \\u_i\\, both
defined on the target grid (`RAS` millimeters). The warped image
returned by
[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md)
samples the source at \\A(r + u_f(r))\\ for every target voxel \\r\\,
and `direction` selects which of the following maps is evaluated:

|  |  |  |  |
|----|----|----|----|
| **call** | **maps** | **composite** | **default output grid** |
| volume, `"forward"` | a source-space image onto the target grid | \\A(r + u_f(r))\\ | target grid |
| volume, `"inverse"` | a target-space image onto the source grid | \\q + u_i(q)\\, \\q = A^{-1} r\\ | source grid |
| points, `"forward"` | source `RAS` to target `RAS` | \\q + u_i(q)\\, \\q = A^{-1} p\\ | (none) |
| points, `"inverse"` | target `RAS` to source `RAS` | \\A(p + u_f(p))\\ | (none) |

Rigid and `affine` registrations behave the same way with \\u \equiv
0\\. With the same `direction`, a volume and a set of points move the
same way (`"forward"`: from source to target space), but the composites
differ: an image is pulled onto its new grid by looking up where each
output voxel comes from, whereas a point is pushed to where it lands, so
the forward point map is the inverse of the forward volume composite.
The inverse field is accurate wherever the forward map is one-to-one
inside the target grid; at the border of an unmasked deformation, where
the field may fold or push points out of the grid, \\u_i\\ has no exact
solution and points mapped there are approximate. The `vox2ras` of
`volume` defaults to the argument, then to the array's `"vox2ras"`
attribute, then to the registration geometry (the source grid for
`"forward"`, the target grid for `"inverse"`); `reference_dim` and
`reference_vox2ras` default to the geometry's target grid for
`"forward"` and source grid for `"inverse"`. A manifest written before
the source dimension was recorded (`SourceDim`) still loads, but then
`reference_dim` must be given for the `"inverse"` direction.
`invert = TRUE` is an error for registration objects: choose the mapping
with `direction`.

## Transform lists

A list in `antsApplyTransforms -t` order. Each element is one of a
\\4\times 4\\ (or \\3\times 4\\) `RAS` matrix mapping fixed to moving
coordinates (such as `$transform` or
[`read_ants_transform`](https://dipterix.org/ravetools/reference/write_ants_transform.md)),
a `(nx, ny, nz, 3)` displacement array in `RAS` millimeters carrying a
`"vox2ras"` attribute (such as `$forward_field` or
[`read_ants_warp`](https://dipterix.org/ravetools/reference/write_ants_warp.md)),
or a file path (`.mat` read with
[`read_ants_transform`](https://dipterix.org/ravetools/reference/write_ants_transform.md),
`.nii` / `.nii.gz` read with
[`read_ants_warp`](https://dipterix.org/ravetools/reference/write_ants_warp.md)).
The composite is \\C = T_n \circ \cdots \circ T_1\\: the *first* element
is applied first to a point, exactly as 'ANTs' does, so transform lists
stored by 'ANTs'-based pipelines work unchanged. Volumes are resampled
as \\out(r) = in(C(r))\\ and points as \\C(p)\\; a registration's
`list(forward_field, transform)` therefore reproduces the `"forward"`
volume warp, and `list(solve(transform), inverse_field)` (or
`list(transform, inverse_field)` with `invert = c(TRUE, FALSE)`) the
`"inverse"` one. `invert` is recycled to the list length and may only be
`TRUE` for matrices; unlike `ANTsPy`, no inversion is ever inferred from
the list layout. A bare matrix or field is accepted as a one-element
list. `reference_dim` and `reference_vox2ras` are required with a list,
and `direction` must not be given.

## Sampling rules

Displacement fields are sampled with `trilinear` interpolation in their
own voxel grid; within the half-voxel margin around their outer nodes
the edge value is used, and beyond that extent they contribute zero
displacement (the rule of the 'ITK' displacement-field transform), so a
point leaving the field is still carried by the remaining stages.
Volumes are sampled with the same rounding and out-of-bounds rules as
[`resample_3d_volume`](https://dipterix.org/ravetools/reference/resample_3d_volume.md)
and
[`apply_transform3d`](https://dipterix.org/ravetools/reference/apply_transform3d.md):
`'nearest'` rounds the continuous voxel coordinate and requires the
index to fall inside the volume, `'trilinear'` requires the coordinate
to lie within the voxel-center box; voxels that do not get `na_fill`.

## See also

[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md),
[`save_registration`](https://dipterix.org/ravetools/reference/save_registration.md),
[`read_ants_transform`](https://dipterix.org/ravetools/reference/write_ants_transform.md),
[`read_ants_warp`](https://dipterix.org/ravetools/reference/write_ants_warp.md),
[`apply_transform3d`](https://dipterix.org/ravetools/reference/apply_transform3d.md)
(superseded, single linear transform only)

## Examples

``` r

# a toy registration: a blob and its shifted copy
nd <- c(24, 24, 24)
vox2ras <- diag(4); vox2ras[1:3, 4] <- -12
blob <- function(cx, cy, cz, s = 4) {
  g <- expand.grid(x = 0:(nd[1]-1), y = 0:(nd[2]-1), z = 0:(nd[3]-1))
  array(exp(-((g$x-cx)^2 + (g$y-cy)^2 + (g$z-cz)^2) / (2*s^2)), nd)
}
target <- blob(12, 12, 12)
source <- blob(14, 11, 12.5)
reg <- register_volume3d(source, target, vox2ras, vox2ras,
                         type = "rigid", metric = "cc", verbose = FALSE)

# warp the source onto the target grid (same as reg$image) ...
warped <- apply_transform3d_volume(source, reg, direction = "forward")
max(abs(warped - reg$image))
#> [1] 0

# ... and a target-space label map back onto the source grid
label <- target > 0.5
label_src <- apply_transform3d_volume(label, reg, direction = "inverse",
                                      interpolation = "nearest")

# points: "forward" maps source RAS to target RAS (like the forward
# volume warp, which brings source-space data into target space)
p_src <- rbind(c(2, -1, 0.5), c(0, 0, 0))
p_tgt <- apply_transform3d_points(p_src, reg, direction = "forward")
apply_transform3d_points(p_tgt, reg, direction = "inverse")   # round trip
#>      [,1]          [,2]          [,3]
#> [1,]    2 -1.000000e+00  5.000000e-01
#> [2,]    0 -1.110223e-16 -1.110223e-16

# the same linear transform as an ANTs-style list (reference grid required)
same <- apply_transform3d_volume(source, list(reg$transform), vox2ras = vox2ras,
                                 reference_dim = nd, reference_vox2ras = vox2ras)
max(abs(same - warped))
#> [1] 0
```
