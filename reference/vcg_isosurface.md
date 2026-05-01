# Create surface mesh from 3D-array

Create surface from 3D-array using marching cubes algorithm

## Usage

``` r
vcg_isosurface(
  volume,
  threshold_lb = 0,
  threshold_ub = NA,
  vox_to_ras = diag(c(-1, -1, 1, 1))
)
```

## Arguments

- volume:

  a volume or a mask volume

- threshold_lb:

  lower-bound threshold for creating the surface; default is `0`

- threshold_ub:

  upper-bound threshold for creating the surface; default is `NA` (no
  upper-bound)

- vox_to_ras:

  a `4x4` `'affine'` transform matrix indicating the 'voxel'-to-world
  transform.

## Value

A triangular mesh of class `'mesh3d'`

## Examples

``` r


if(is_not_cran()) {

library(ravetools)
data("left_hippocampus_mask")

mesh <- vcg_isosurface(left_hippocampus_mask)


rgl_view({

  rgl_call("mfrow3d", 1, 2)

  rgl_call("title3d", "Direct ISOSurface")
  rgl_call("shade3d", mesh, col = 2)

  rgl_call("next3d")
  rgl_call("title3d", "ISOSurface + Implicit Smooth")

  rgl_call("shade3d",
           vcg_smooth_implicit(mesh, degree = 2),
           col = 3)
})

}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’
```
