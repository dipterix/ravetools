# Sub-divide (up-sample) a triangular mesh

Up-sample a triangular mesh by adding a vertex at each edge or face
center.

## Usage

``` r
vcg_subdivision(mesh, method = c("edge", "barycenter"))
```

## Arguments

- mesh:

  triangular mesh stored as object of class 'mesh3d'.

- method:

  either `'edge'` (default) to add new mid-point vertices to edge, or
  `'barycenter'` to add new vertices at face `'Bary'` centers.

## Value

An object of class "mesh3d"

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

mesh <- plane_geometry()

# default
mesh_edge <- vcg_subdivision(mesh, "edge")

# barycenter
mesh_face <- vcg_subdivision(mesh, "barycenter")

if(is_not_cran()) {

  rgl_view({
    rgl_call("wire3d", mesh, col = 1)
    rgl_call("wire3d", mesh_edge, col = 2)
    rgl_call("wire3d", mesh_face, col = 3)
  })


}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’


```
