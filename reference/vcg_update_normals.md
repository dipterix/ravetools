# Update vertex normal

Update vertex normal

## Usage

``` r
vcg_update_normals(
  mesh,
  weight = c("area", "angle"),
  pointcloud = c(10, 0),
  verbose = FALSE
)
```

## Arguments

- mesh:

  triangular mesh or a point-cloud (matrix of 3 columns)

- weight:

  method to compute per-vertex normal vectors: `"area"` weighted average
  of surrounding face normal, or `"angle"` weighted vertex normal
  vectors.

- pointcloud:

  integer vector of length 2: containing optional parameters for normal
  calculation of point clouds; the first entry specifies the number of
  neighboring points to consider; the second entry specifies the amount
  of smoothing iterations to be performed.

- verbose:

  whether to verbose the progress

## Value

A `'mesh3d'` object with normal vectors.

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

if(is_not_cran()) {

# Prepare mesh with no normal
data("left_hippocampus_mask")
mesh <- vcg_isosurface(left_hippocampus_mask)
mesh$normals <- NULL

# Start: examples
new_mesh <- vcg_update_normals(mesh, weight = "angle",
                               pointcloud = c(10, 10))

rgl_view({
  rgl_call("mfrow3d", 1, 2)
  rgl_call("shade3d", mesh, col = 2)

  rgl_call("next3d")
  rgl_call("shade3d", new_mesh, col = 2)
})
}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’

```
