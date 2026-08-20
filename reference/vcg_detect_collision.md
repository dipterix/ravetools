# Detect collisions between two geometries

Reports, for every element of `y`, whether it comes within `radius` of
`x`, together with the exact minimum distance. Each side may be a point
cloud, a set of connected line segments (for example diffusion
streamlines or an electrode shaft), or a triangular mesh, so the
function covers all nine pairings with one call.

The test is exact: it produces neither false negatives nor false
positives, and `distance` is the true minimum distance rather than a
vertex-sampled approximation.

## Usage

``` r
vcg_detect_collision(
  x,
  y,
  mode_x = c("auto", "points", "segments", "mesh"),
  mode_y = c("auto", "points", "segments", "mesh"),
  radius = 0,
  early_stop = FALSE,
  include_interior = FALSE
)
```

## Arguments

- x:

  the geometry that gets indexed; a matrix with `n` rows and 2 or 3
  columns, or a surface (see
  [`ensure_mesh3d`](https://dipterix.org/ravetools/reference/ensure_mesh3d.md))

- y:

  the geometry that gets queried; results are aligned to this side. Same
  accepted types as `x`

- mode_x, mode_y:

  how to interpret `x` and `y`: `'points'` treats each row as an
  isolated point, `'segments'` chains consecutive rows together
  (`p1-p2`, `p2-p3`, ...), and `'mesh'` uses the triangular faces. The
  default `'auto'` resolves to `'mesh'` for a surface that has faces and
  `'points'` otherwise, so a matrix of connected line segments must ask
  for `'segments'` explicitly

- radius:

  distance tolerance; a collision is reported when the minimum distance
  is at most `radius`. Defaults to `0`, meaning literal contact; a
  positive tolerance is usually what is wanted

- early_stop:

  whether to stop as soon as a collision is found, scanning
  independently within each group of `y` delimited by `NA` rows. A
  matrix carrying no `NA` rows is a single group, so a point cloud stops
  after one hit overall rather than per row; likewise `mode_y = 'mesh'`
  stops at the first colliding face. Skipped elements report `hit = NA`,
  meaning *not evaluated* rather than *no collision*, and the collision
  reported is the first in row order, not the closest. Use `FALSE` (the
  default) whenever the per-element profile or a minimum distance is
  wanted

- include_interior:

  whether geometry lying strictly inside a closed `x` surface counts as
  a collision even when it never comes within `radius` of the surface
  itself; requires `mode_x = 'mesh'`. Default is `FALSE`

## Value

A list of three vectors, each aligned one-to-one with the rows of `y`
(or with the faces of `y` when `mode_y = 'mesh'`):

- `hit`:

  logical. `NA` marks an element with nothing to test: a separator row,
  or the final vertex of each group when `mode_y = 'segments'`. When
  `early_stop = TRUE`, `NA` also marks elements the scan never reached

- `distance`:

  the exact minimum distance where `hit` is `TRUE`, and `NA` everywhere
  else. Distances beyond `radius` are never computed, so a `FALSE`
  element carries no distance

- `index`:

  1-based index of the closest element of `x` within `radius`: a row
  when `mode_x` is `'points'`, the row where the segment starts when
  `'segments'`, or a face column when `'mesh'`. `NA` where there is no
  collision

## Details

`x` is indexed once into a uniform spatial grid and `y` is streamed
against it, so put the larger or repeatedly-reused geometry in `x` and
the geometry you want per-element answers about in `y`.

Several chains share one matrix and are delimited by rows of `NA`; any
row that is not fully finite is treated as a separator. This lets a
whole bundle of streamlines be passed as a single matrix:


      rbind(
        c(0, 0, 0), c(1, 0, 0), c(2, 0, 0),   # first streamline
        c(NA, NA, NA),                        # separator
        c(0, 5, 0), c(1, 5, 0)                # second streamline
      )

Those separators also define the unit that `early_stop` works on: it
stops once per group, so a bundle of streamlines yields at most one hit
per streamline.

When `mode_y = 'segments'`, row `i` of the result describes the segment
running from row `i` to row `i + 1`.

Both geometries must already share one coordinate space. Surface
coordinates, volume `IJK` indices, and scanner coordinates are all
different spaces, and no transform is applied here.

`include_interior` relies on ray casting and therefore needs `x` to be
watertight with coherently oriented faces; see
[`vcg_fix_defects`](https://dipterix.org/ravetools/reference/vcg_fix_defects.md)
to repair a surface that is not.

Multi-threading follows
[`ravetools_threads`](https://dipterix.org/ravetools/reference/parallel-options.md).
The interior test always runs single-threaded.

## Coercing Surface Inputs

The surface objects are converted to `'mesh3d'` object before applying
further calculations.

When `surface` is a surface ieegio object, the returned `mesh3d$vb`
contains vertices that have been left-multiplied by
`surface$geometry$transforms[[1]]` (the first transform stored in the
geometry, typically the `ScannerAnat` or voxel-to-world transform).

**Breaking change:** Earlier versions (before 0.2.6) of ravetools
returned the raw `surface$geometry$vertices` without applying any
transform, so downstream code often multiplied by
`surface$geometry$transforms[[1]]` (or an equivalent) manually before
working in world space. Such code will now *double* apply the transform
and produce incorrect coordinates. If you previously applied a transform
from `surface$geometry$transforms` by hand after calling a ravetools
mesh function on an `'ieegio_surface'`, remove that manual step.

Surfaces with an empty or missing `geometry$transforms` list (for
example, surfaces produced by ieegio's `volume_to_surface`, which stores
an identity transform) are unaffected.

If `geometry$transforms` contains multiple transforms targeting
different coordinate spaces, only the first one is used. Callers that
need a specific target space should select and apply that transform
themselves before calling ravetools mesh functions.

## Examples

``` r

library(ravetools)

# A spherical region of interest
roi <- vcg_sphere()

# Two streamlines in one matrix, separated by an NA row: the first passes
# through the sphere, the second stays well clear of it
streamlines <- rbind(
  cbind(seq(-3, 3, by = 0.5), 0, 0),
  c(NA, NA, NA),
  cbind(seq(-3, 3, by = 0.5), 5, 0)
)

result <- vcg_detect_collision(
  x = roi,
  y = streamlines,
  mode_y = "segments",
  radius = 0.1
)

# Which segments touch the sphere, and how close do they get?
data.frame(
  hit = result$hit,
  distance = round(result$distance, 4)
)
#>      hit distance
#> 1  FALSE       NA
#> 2  FALSE       NA
#> 3  FALSE       NA
#> 4   TRUE        0
#> 5   TRUE        0
#> 6  FALSE       NA
#> 7  FALSE       NA
#> 8   TRUE        0
#> 9   TRUE        0
#> 10 FALSE       NA
#> 11 FALSE       NA
#> 12 FALSE       NA
#> 13    NA       NA
#> 14    NA       NA
#> 15 FALSE       NA
#> 16 FALSE       NA
#> 17 FALSE       NA
#> 18 FALSE       NA
#> 19 FALSE       NA
#> 20 FALSE       NA
#> 21 FALSE       NA
#> 22 FALSE       NA
#> 23 FALSE       NA
#> 24 FALSE       NA
#> 25 FALSE       NA
#> 26 FALSE       NA
#> 27    NA       NA

# A point cloud against the same region, asking only for proximity
points <- rbind(c(0, 0, 0), c(1.05, 0, 0), c(10, 10, 10))
vcg_detect_collision(roi, points, radius = 0.1)$hit
#> [1] FALSE  TRUE FALSE

# The centre of the sphere is far from its surface, so it only counts as a
# collision when the interior is included
vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1)$hit
#> [1] FALSE
vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1,
                     include_interior = TRUE)$hit
#> [1] TRUE
```
