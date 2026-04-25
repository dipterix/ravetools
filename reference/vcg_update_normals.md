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
