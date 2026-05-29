# Coerce a surface object to a `'mesh3d'` mesh

Internal helper used throughout ravetools to accept a variety of surface
representations and return a list of class `'mesh3d'` (the format used
by the `rgl` package). Most mesh-consuming functions in this package
call `ensure_mesh3d` on their surface arguments, so the coercion rules
described here apply to all of them.

## Usage

``` r
ensure_mesh3d(surface)
```

## Arguments

- surface:

  a surface object. One of the following:

  `'mesh3d'`

  :   returned unchanged.

  `'fs.surface'`

  :   a surface in the format produced by
      [`freesurferformats::read.fs.surface`](https://rdrr.io/pkg/freesurferformats/man/read.fs.surface.html).
      Vertices and faces are copied into a `'mesh3d'` list; zero-indexed
      faces are bumped by 1.

  `'ieegio_surface'`

  :   a surface object produced by ieegio (e.g. from `read_surface` or
      `volume_to_surface`). See the section below for important behavior
      changes.

  a bare `list`

  :   must contain a numeric `vb` matrix with 3 or 4 rows; the list is
      reclassified as `'mesh3d'` with no other modification.

  Any other input triggers an error.

## Value

An object of class `'mesh3d'` with at least the `vb` (vertex) component
and, when face information is available, an `it` (triangle index)
component.

## Coercing `ieegio_surface` inputs

When `surface` is an `'ieegio_surface'` object, the returned `mesh3d$vb`
contains vertices that have been left-multiplied by
`surface$geometry$transforms[[1]]` (the first transform stored in the
geometry, typically the `ScannerAnat` or voxel-to-world transform).

**Breaking change:** Earlier versions of ravetools returned the raw
`surface$geometry$vertices` without applying any transform, so
downstream code often multiplied by `surface$geometry$transforms[[1]]`
(or an equivalent) manually before working in world space. Such code
will now *double* apply the transform and produce incorrect coordinates.
If you previously applied a transform from `surface$geometry$transforms`
by hand after calling a ravetools mesh function on an
`'ieegio_surface'`, remove that manual step.

Surfaces with an empty or missing `geometry$transforms` list (for
example, surfaces produced by ieegio's `volume_to_surface`, which stores
an identity transform) are unaffected.

If `geometry$transforms` contains multiple transforms targeting
different coordinate spaces, only the first one is used. Callers that
need a specific target space should select and apply that transform
themselves before calling ravetools mesh functions.

## Examples

``` r

# mesh3d input is returned unchanged in shape
sphere <- vcg_sphere()
m <- ensure_mesh3d(sphere)
identical(m$vb, sphere$vb)
#> [1] TRUE

# A bare list with a `vb` slot is reclassified to mesh3d
bare <- list(vb = rbind(matrix(rnorm(30), nrow = 3), 1))
m <- ensure_mesh3d(bare)
inherits(m, "mesh3d")
#> [1] TRUE
```
