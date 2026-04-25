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
