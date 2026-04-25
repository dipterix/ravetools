# Generate 3D mesh surface from volume data

This function is soft-deprecated. Please use
[`vcg_mesh_volume`](https://dipterix.org/ravetools/reference/vcg_mesh_volume.md),
[`vcg_uniform_remesh`](https://dipterix.org/ravetools/reference/vcg_uniform_remesh.md),
and
[`vcg_smooth_explicit`](https://dipterix.org/ravetools/reference/vcg_smooth.md)
or
[`vcg_smooth_implicit`](https://dipterix.org/ravetools/reference/vcg_smooth.md).

## Usage

``` r
mesh_from_volume(
  volume,
  output_format = c("rgl", "freesurfer"),
  IJK2RAS = NULL,
  threshold = 0,
  verbose = TRUE,
  remesh = TRUE,
  remesh_voxel_size = 1,
  remesh_multisample = TRUE,
  remesh_automerge = TRUE,
  smooth = FALSE,
  smooth_lambda = 10,
  smooth_delta = 20,
  smooth_method = "surfPreserveLaplace"
)
```

## Arguments

- volume:

  3-dimensional volume array

- output_format:

  resulting data format, choices are `'rgl'` and `'freesurfer'`

- IJK2RAS:

  volume 'IJK' (zero-indexed coordinate index) to `'tkrRAS'` transform,
  default is automatically determined

- threshold:

  threshold used to create volume mask; the surface will be created to
  fit the mask boundaries

- verbose:

  whether to verbose the progress

- remesh:

  whether to re-sample the mesh using
  [`vcg_uniform_remesh`](https://dipterix.org/ravetools/reference/vcg_uniform_remesh.md)

- remesh_voxel_size, remesh_multisample, remesh_automerge:

  see arguments in
  [`vcg_uniform_remesh`](https://dipterix.org/ravetools/reference/vcg_uniform_remesh.md)

- smooth:

  whether to smooth the mesh via
  [`vcg_smooth_explicit`](https://dipterix.org/ravetools/reference/vcg_smooth.md)

- smooth_lambda, smooth_delta, smooth_method:

  see
  [`vcg_smooth_explicit`](https://dipterix.org/ravetools/reference/vcg_smooth.md)

## Value

A `'mesh3d'` surface if `output_format` is 'rgl', or `'fs.surface'`
surface otherwise.

## Examples

``` r

volume <- array(0, dim = c(8,8,8))
volume[4:5, 4:5, 4:5] <- 1

graphics::image(x = volume[4,,])


# you can use rgl::wire3d(mesh) to visualize the mesh
mesh <- mesh_from_volume(volume, verbose = FALSE)

```
