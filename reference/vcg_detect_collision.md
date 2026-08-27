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
  test_level = c("element", "unit", "whole"),
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

- test_level:

  how exhaustively `y` is scanned. Answers are always reported per
  *unit* of `y`: a point when `mode_y = 'points'`, a chain when
  `'segments'`, a face when `'mesh'`. `'element'` (the default) measures
  every element of every unit and reports the closest one; `'unit'`
  stops at the first hit inside each unit and reports that one;
  `'whole'` stops at the first hit anywhere and answers only whether the
  two geometries collide. Because a point and a face each *are* their
  own unit, `'element'` and `'unit'` differ only when
  `mode_y = 'segments'`

- include_interior:

  whether geometry lying strictly inside a closed `x` surface counts as
  a collision even when it never comes within `radius` of the surface
  itself; requires `mode_x = 'mesh'`. Default is `FALSE`

## Value

A list. `collide` is always present, so a plain yes-or-no question can
be asked at any `test_level`:

- `collide`:

  a length-one logical: whether `x` and `y` collide at all. `FALSE` for
  empty input

- `hit_unit`:

  logical, one entry per unit of `y`, omitted when
  `test_level = 'whole'`. `NA` marks a unit with nothing to test, and is
  deliberately not `FALSE`; see *Why `hit_unit` reports `NA`* below

- `representation`:

  a data frame with one row per *hit* unit, so `nrow` equals
  `sum(hit_unit, na.rm = TRUE)` and a collision-free query gives zero
  rows. When `test_level = 'whole'` it holds at most one row

- `summary`:

  what the query actually was, so a stored result can be read back
  without the call that made it

`summary` holds `radius`, `test_level`, `include_interior`, and one
entry per side. Each side reports its `mode` *as resolved* (never
`'auto'`), the matching `unit_type` of `'point'`, `'line'`, or `'face'`,
and `n_units`: one per row for `'points'` with separator rows included,
one per chain for `'segments'`, one per face for `'mesh'`. So
`summary$y$n_units` always equals `length(hit_unit)`, and remains the
only record of the size of `y` when `test_level = 'whole'` drops that
vector.

The columns of `representation` name the collision on both sides:

- `unit`:

  1-based unit of `y`: the row itself when `mode_y = 'points'`, the
  position of the chain when `'segments'`, the face column when `'mesh'`

- `index`:

  1-based element *within* that unit: always `1` for a point, the
  position of the segment counted from the start of its own chain for
  `'segments'`, and the face vertex nearest `x` (`1`, `2`, or `3`) for
  `'mesh'`

- `distance`:

  the exact distance from that element to `x`. Distances beyond `radius`
  are never computed

- `x_index`:

  1-based element of `x` that was hit: a row when `mode_x` is
  `'points'`, the row where the segment starts when `'segments'`, or a
  face column when `'mesh'`. `NA` for a hit found by `include_interior`
  rather than by proximity

## Details

`x` is indexed once into a uniform spatial grid and `y` is streamed
against it, so put the larger or repeatedly-reused geometry in `x` and
the geometry you want per-unit answers about in `y`.

Several chains share one matrix and are delimited by rows of `NA`; any
row that is not fully finite is treated as a separator. This lets a
whole bundle of streamlines be passed as a single matrix:


      rbind(
        c(0, 0, 0), c(1, 0, 0), c(2, 0, 0),   # first streamline
        c(NA, NA, NA),                        # separator
        c(0, 5, 0), c(1, 5, 0)                # second streamline
      )

Those separators define the units that results are reported against:
with `mode_y = 'segments'` the matrix above is two units, so `hit_unit`
has two entries whatever the `test_level`. With `mode_y = 'points'`
every row is its own unit, separators included, so `hit_unit` stays
aligned with the rows of `y` and reports `NA` at each separator.

Separators split the matrix into slots, and *every* slot is a chain,
empty ones included: two adjacent separators enclose an empty chain,
which reports `NA`. The two ends are deliberately not symmetric. A
separator in the first row opens an empty leading chain, but a separator
in the last row closes the chain before it rather than opening a
trailing empty one:

|             |                 |
|-------------|-----------------|
| **matrix**  | **chains**      |
| `A`         | 1: `[A]`        |
| `A NA B`    | 2: `[A] [B]`    |
| `A NA NA B` | 3: `[A] [] [B]` |
| `NA A`      | 2: `[] [A]`     |
| `A NA`      | 1: `[A]`        |

That rule exists so a bundle assembled by appending a separator to every
chain gives exactly one unit per chain, whichever of them happen to be
empty:


      do.call("rbind", lapply(tracts, function(line) rbind(line, NA)))

Since `mode_x` plays no part in that grouping, the `index` of a segment
is always counted from the start of its own chain rather than from the
top of the matrix. The `unit` column gives the chain.

## Why `hit_unit` reports `NA`

`NA` means *there was no geometry here to test*. `FALSE` means *geometry
was here and it missed*. Three situations produce the first:

|              |                                                   |
|--------------|---------------------------------------------------|
| **`mode_y`** | **the unit reports `NA` when**                    |
| `'points'`   | the row is a separator                            |
| `'segments'` | the chain is empty (two adjacent separators)      |
| `'segments'` | the chain holds one vertex, so carries no segment |

`mode_y = 'mesh'` therefore never reports `NA`: a face is always
testable. (The scan does skip faces flagged as deleted, but a surface is
rebuilt from its vertex and face matrices on every call, so no face
arrives flagged.)

An empty chain cannot collide with anything, so reporting `FALSE` would
not be wrong in the set-theoretic sense. It is avoided for two reasons.
In `'points'` mode `hit_unit` is aligned to the rows of `y` so that it
can be indexed straight back into the input matrix, which means
separator rows occupy slots; calling those `FALSE` asserts something
about a point that does not exist. And `sum(is.na(hit_unit))` counts the
degenerate units that were passed in – empty or single-vertex
streamlines in a tract file, deleted faces in a surface – which is
usually a data-quality signal worth seeing rather than silently folding
into "missed". Once flattened to `FALSE` that count cannot be recovered.

Most uses never have to deal with it. `which` already drops `NA`, so
`which(hit_unit)` gives the colliding units directly, and
`representation$unit` is that same vector without any coercion. `NA`
only intrudes on `sum` / `any` / `if` without `na.rm`, and on logical
subsetting such as `tracts[hit_unit]`, which yields `NULL` entries.
Where a plain mask really is wanted:


      mask <- res$hit_unit %in% TRUE     # NA becomes FALSE

`test_level = 'whole'` reports *a* colliding unit, not necessarily the
first one in order: the scan stops as soon as any thread finds a
collision. `collide` itself is deterministic, as are the two lower
levels, which measure each unit in row order.

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

# Which streamlines touch the sphere, and where does each come closest?
result$hit_unit
#> [1]  TRUE FALSE
result$representation
#>   unit index distance x_index
#> 1    1     4        0     494

# What was actually asked, including what `mode = "auto"` resolved to
result$summary
#> $x
#> $x$mode
#> [1] "mesh"
#> 
#> $x$unit_type
#> [1] "face"
#> 
#> $x$n_units
#> [1] 1280
#> 
#> 
#> $y
#> $y$mode
#> [1] "segments"
#> 
#> $y$unit_type
#> [1] "line"
#> 
#> $y$n_units
#> [1] 2
#> 
#> 
#> $radius
#> [1] 0.1
#> 
#> $test_level
#> [1] "element"
#> 
#> $include_interior
#> [1] FALSE
#> 

# A bundle assembled by appending a separator to every streamline gives one
# unit per streamline, even where a streamline is empty: here the middle one
# contributes only its own separator, and reports NA rather than vanishing
tracts <- list(
  cbind(seq(-3, 3, by = 0.5), 0, 0),
  matrix(numeric(0), ncol = 3),
  cbind(seq(-3, 3, by = 0.5), 5, 0)
)
bundle <- do.call("rbind", lapply(tracts, function(line) rbind(line, NA)))
vcg_detect_collision(roi, bundle, mode_y = "segments",
                     radius = 0.1)$hit_unit
#> [1]  TRUE    NA FALSE

# The same question, giving up on a streamline once it is known to collide.
# `hit_unit` is unchanged; `index` names the first colliding segment instead
# of the closest one, which here happens to be the same segment
vcg_detect_collision(roi, streamlines, mode_y = "segments", radius = 0.1,
                     test_level = "unit")$representation
#>   unit index distance x_index
#> 1    1     4        0     494

# Does the bundle touch the sphere at all?
vcg_detect_collision(roi, streamlines, mode_y = "segments", radius = 0.1,
                     test_level = "whole")$collide
#> [1] TRUE

# A point cloud against the same region, asking only for proximity
points <- rbind(c(0, 0, 0), c(1.05, 0, 0), c(10, 10, 10))
vcg_detect_collision(roi, points, radius = 0.1)$hit_unit
#> [1] FALSE  TRUE FALSE

# The centre of the sphere is far from its surface, so it only counts as a
# collision when the interior is included
vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1)$collide
#> [1] FALSE
vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1,
                     include_interior = TRUE)$collide
#> [1] TRUE
```
