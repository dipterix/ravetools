# Project plane to a surface

Project a two-dimensional plane (such as `'ECoG'` grid) to a
three-dimensional surface while preserving the order

## Usage

``` r
project_plane(
  target,
  width,
  height,
  shape,
  initial_positions,
  translate_first = TRUE,
  diagnostic = FALSE,
  n_iters = 5
)
```

## Arguments

- target:

  target surface to be projected to, must be object that can be
  converted to `'mesh3d'` (`'rgl'` surface), for example, `'fs.surface'`
  (from `'freesurferformat'` package) or `'ieegio_surface'` from ieegio
  package; see
  [`ensure_mesh3d`](https://dipterix.org/ravetools/reference/ensure_mesh3d.md)
  for the coercion rules.

- width, height:

  width and height of the plane in world space (for `'ECoG'` grid, the
  unit is millimeter)

- shape:

  vector of two integers: the first element is the number of vertices
  (or electrode contacts) along `'width'` direction; the second element
  is the number of vertices along `'height'` direction. The total number
  of vertices of the plane will be `prod(shape)`. Notice

- initial_positions:

  a `shape[[1]] x shape[[2]] x 3` array or a `n x 3` matrix, where `n`
  is `prod(shape)`, the number of vertices indicating the initial vertex
  positions of the plane

- translate_first:

  whether to translate the plane first if the plane center is far from
  the surface; default is `FALSE`; set to `TRUE` for a warm start

- diagnostic:

  whether to plot diagnostic figures showing the morphing progress.

- n_iters:

  number of iterations; default is five

## Value

The projected vertex locations, same order as `initial_positions`.

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


# Construct target surface

sphere <- vcg_sphere()
target <- structure(
  class = "mesh3d",
  list(
    vb = cbind(
      sphere$vb[1:3, ] - c(0.8, 0, 0),
      sphere$vb[1:3, ] + c(0.8, 0, 0)
    ),
    it = cbind(
      sphere$it[1:3, ],
      sphere$it[1:3, ] + ncol(sphere$vb)
    )
  )
)
n_surfverts <- ncol(target$vb)

plane <- plane_geometry(width = 3, height = 3, shape = c(30, 30))
plane$vb <- plane$vb[1:3, , drop = FALSE] + c(0, 0, 2)
n_contacts <- ncol(plane$vb)

# First plot
x <- t(cbind(target$vb, plane$vb))
colnames(x) <- c('x', 'y', 'z')
graphics::pairs(
  x = x, asp = 1,
  col = c(
    rep("black", n_surfverts),
    rep("green", n_contacts)
  ),
  pch = c(
    rep(46, n_surfverts),
    rep(20, n_contacts)
  )
)


projected <- project_plane(
  target = target, width = 3, height = 3, shape = c(30, 30),
  initial_positions = t(plane$vb),
  translate_first = TRUE, diagnostic = FALSE
)

y <- rbind(x, projected)
graphics::pairs(
  x = y, asp = 1,
  col = c(
    rep("black", ncol(target$vb)),
    rep("green", n_contacts),
    rep("red", n_contacts)
  ),
  pch = c(
    rep(46, n_surfverts),
    rep(1, n_contacts),
    rep(20, n_contacts)
  )
)


```
