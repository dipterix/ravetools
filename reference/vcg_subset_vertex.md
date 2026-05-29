# Subset mesh by vertex

Subset mesh by vertex

## Usage

``` r
vcg_subset_vertex(x, selector)
```

## Arguments

- x:

  surface mesh

- selector:

  logical vector (must not contain NA), and length must be consistent
  with the number of vertices in `x`: which nodes are to be kept

## Value

A triangular mesh of class `'mesh3d'`, a subset of `x`

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

sphere <- vcg_sphere()

nv <- ncol(sphere$vb)

selector <- seq_len(nv) > (nv / 2)

sub <- vcg_subset_vertex(sphere, selector)

if(is_not_cran()) {
  rgl_view({

    # subset sphere will be displayed in red
    rgl_call("shade3d", sub, col = 'red')

    # Original sphere will be displayed as wireframe
    rgl_call("wire3d", sphere, col = (2 - selector))

  })
}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’

```
