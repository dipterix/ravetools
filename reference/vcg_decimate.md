# Simplify a triangular mesh by `quadric` edge collapse

Reduces the number of faces of a triangular mesh by repeatedly
collapsing the edge whose removal changes the surface least, measured by
the `quadric` error metric of `Garland` and `Heckbert` (1997): the
squared distance from the merged vertex to the planes of the faces
around it. Flat regions are simplified first and curved ones last, so
the shape is kept while the vertex count drops. Useful before
[`vcg_smooth_implicit`](https://dipterix.org/ravetools/reference/vcg_smooth.md)
on very large meshes, such as surfaces extracted from whole-brain
volumes.

## Usage

``` r
vcg_decimate(
  mesh,
  ratio = 0.5,
  target_faces = NULL,
  preserve_topology = TRUE,
  preserve_boundary = TRUE,
  normal_check = TRUE,
  quality_threshold = 0.3,
  verbose = FALSE
)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`

- ratio:

  fraction of the faces to keep, greater than 0 and at most 1; ignored
  when `target_faces` is given

- target_faces:

  number of faces to keep; default `NULL` uses `ratio`

- preserve_topology:

  whether to forbid collapses that change the topology (join or split
  connected parts, open or close handles); default is `TRUE`

- preserve_boundary:

  whether to keep the boundary edges of an open mesh in place; default
  is `TRUE`

- normal_check:

  whether to forbid collapses that flip a face; default is `TRUE`

- quality_threshold:

  collapses that would create faces of lower quality than this, from 0
  (degenerate) to 1 (equilateral), are penalized; default is `0.3`

- verbose:

  whether to print the vertex and face counts

## Value

A `'mesh3d'` object with `vb`, `it`, and `normals`; vertices are
re-indexed. A mesh that already has no more faces than requested is
returned unchanged. The decimation can stop short of the target when the
constraints above forbid every remaining collapse.

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

## References

`Garland` M, `Heckbert` PS (1997). Surface simplification using
`quadric` error metrics. In *Proceedings of the annual conference on
computer graphics and interactive techniques*, 209-216.

## Examples

``` r

sphere <- vcg_sphere(sub_division = 4L)
ncol(sphere$it)
#> [1] 5120

simplified <- vcg_decimate(sphere, ratio = 0.1)
ncol(simplified$it)
#> [1] 512
```
