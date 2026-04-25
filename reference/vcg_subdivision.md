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
