# Implicitly smooth a triangular mesh

Applies smoothing algorithms on a triangular mesh.

## Usage

``` r
vcg_smooth_implicit(
  mesh,
  lambda = 0.2,
  use_mass_matrix = TRUE,
  fix_border = FALSE,
  use_cot_weight = FALSE,
  degree = 1L,
  laplacian_weight = 1
)

vcg_smooth_explicit(
  mesh,
  type = c("taubin", "laplace", "HClaplace", "fujiLaplace", "angWeight",
    "surfPreserveLaplace"),
  iteration = 10,
  lambda = 0.5,
  mu = -0.53,
  delta = 0.1
)
```

## Arguments

- mesh:

  triangular mesh stored as object of class 'mesh3d'.

- lambda:

  In `vcg_smooth_implicit`, the amount of smoothness, useful only if
  `use_mass_matrix` is `TRUE`; default is `0.2`. In
  `vcg_smooth_explicit`, parameter for `'taubin'` smoothing.

- use_mass_matrix:

  logical: whether to use mass matrix to keep the mesh close to its
  original position (weighted per area distributed on vertices); default
  is `TRUE`

- fix_border:

  logical: whether to fix the border vertices of the mesh; default is
  `FALSE`

- use_cot_weight:

  logical: whether to use cotangent weight; default is `FALSE` (using
  uniform 'Laplacian')

- degree:

  integer: degrees of 'Laplacian'; default is `1`

- laplacian_weight:

  numeric: weight when `use_cot_weight` is `FALSE`; default is `1.0`

- type:

  method name of explicit smooth, choices are `'taubin'`, `'laplace'`,
  `'HClaplace'`, `'fujiLaplace'`, `'angWeight'`,
  `'surfPreserveLaplace'`.

- iteration:

  number of iterations

- mu:

  parameter for `'taubin'` explicit smoothing.

- delta:

  parameter for scale-dependent 'Laplacian' smoothing or maximum allowed
  angle (in 'Radian') for deviation between surface preserving
  'Laplacian'.

## Value

An object of class "mesh3d" with:

- `vb`:

  vertex coordinates

- `normals`:

  vertex normal vectors

- `it`:

  triangular face index

## Examples

``` r
if(is_not_cran()) {

# Prepare mesh with no normals
data("left_hippocampus_mask")

# Grow 2mm on each direction to fill holes
volume <- grow_volume(left_hippocampus_mask, 2)

# Initial mesh
mesh <- vcg_isosurface(volume)

# Start: examples
rgl_view({
  rgl_call("mfrow3d", 2, 4)
  rgl_call("title3d", "Naive ISOSurface")
  rgl_call("shade3d", mesh, col = 2)

  rgl_call("next3d")
  rgl_call("title3d", "Implicit Smooth")
  rgl_call("shade3d", col = 2,
           x = vcg_smooth_implicit(mesh, degree = 2))

  rgl_call("next3d")
  rgl_call("title3d", "Explicit Smooth - taubin")
  rgl_call("shade3d", col = 2,
           x = vcg_smooth_explicit(mesh, "taubin"))

  rgl_call("next3d")
  rgl_call("title3d", "Explicit Smooth - laplace")
  rgl_call("shade3d", col = 2,
           x = vcg_smooth_explicit(mesh, "laplace"))

  rgl_call("next3d")
  rgl_call("title3d", "Explicit Smooth - angWeight")
  rgl_call("shade3d", col = 2,
           x = vcg_smooth_explicit(mesh, "angWeight"))

  rgl_call("next3d")
  rgl_call("title3d", "Explicit Smooth - HClaplace")
  rgl_call("shade3d", col = 2,
           x = vcg_smooth_explicit(mesh, "HClaplace"))

  rgl_call("next3d")
  rgl_call("title3d", "Explicit Smooth - fujiLaplace")
  rgl_call("shade3d", col = 2,
           x = vcg_smooth_explicit(mesh, "fujiLaplace"))

  rgl_call("next3d")
  rgl_call("title3d", "Explicit Smooth - surfPreserveLaplace")
  rgl_call("shade3d", col = 2,
           x = vcg_smooth_explicit(mesh, "surfPreserveLaplace"))
})

}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’
```
