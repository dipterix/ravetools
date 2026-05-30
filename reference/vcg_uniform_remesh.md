# Sample a surface mesh uniformly

Sample a surface mesh uniformly

## Usage

``` r
vcg_uniform_remesh(
  x,
  voxel_size = NULL,
  offset = 0,
  discretize = FALSE,
  multi_sample = FALSE,
  absolute_distance = FALSE,
  merge_clost = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  surface

- voxel_size:

  'voxel' size for space 'discretization'

- offset:

  offset position shift of the new surface from the input

- discretize:

  whether to use step function(`TRUE`) instead of linear interpolation
  (`FALSE`) to calculate the position of the intersected edge of the
  marching cube; default is `FALSE`

- multi_sample:

  whether to calculate multiple samples for more accurate results (at
  the expense of more computing time) to remove artifacts; default is
  `FALSE`

- absolute_distance:

  whether an unsigned distance field should be computed. When set to
  `TRUE`, non-zero offsets is to be set, and double-surfaces will be
  built around the original surface, like a sandwich.

- merge_clost:

  whether to merge close vertices; default is `TRUE`

- verbose:

  whether to verbose the progress; default is `TRUE`

## Value

A triangular mesh of class `'mesh3d'`

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

sphere <- vcg_sphere()
mesh <- vcg_uniform_remesh(sphere, voxel_size = 0.45)
#> Resampling mesh using a volume of 5 x 5 x 5
#>   VoxelSize is 0.450000, offset is 0.000000
#>   Mesh Box is 2.000000 2.000000 2.000000

if(is_not_cran()) {

rgl_view({

  rgl_call("mfrow3d", 1, 2)

  rgl_call("title3d", "Input")
  rgl_call("wire3d", sphere, col = 2)
  rgl_call("next3d")

  rgl_call("title3d", "Re-meshed to 0.1mm edge distance")
  rgl_call("wire3d", mesh, col = 3)
})

}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’
```
