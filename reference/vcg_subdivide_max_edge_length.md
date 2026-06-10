# Selectively subdivide mesh edges that exceed a length threshold

Up-sample a triangular mesh by iteratively splitting only edges longer
than `max_edge_len`. Each long edge is split at its midpoint; the new
vertex is connected to the opposite corner of every adjacent face.
Iteration stops when no edge exceeds the threshold or `max_iter` passes
are exhausted.

This is far cheaper than
[`vcg_subdivision`](https://dipterix.org/ravetools/reference/vcg_subdivision.md)
when most edges are already short and only a small fraction need
splitting.

## Usage

``` r
vcg_subdivide_max_edge_length(mesh, max_edge_len, max_iter = NULL)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`.

- max_edge_len:

  maximum allowed edge length (same units as mesh coordinates).

- max_iter:

  maximum number of refinement passes. When `NULL` (default), derived
  automatically from the current maximum edge length:
  `ceiling(log2(current_max / max_edge_len)) + 1L`, the minimum number
  of bisections needed in the worst case. The loop also terminates early
  once no edge exceeds the threshold.

## Value

An object of class `"mesh3d"` with all edges at most `max_edge_len` long
(provided `max_iter` was sufficient).

## Note

The mesh must be manifold. Run
[`vcg_fix_defects`](https://dipterix.org/ravetools/reference/vcg_fix_defects.md)
first if the mesh has boundary edges or non-manifold vertices.

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

## See also

[`vcg_max_edge_length`](https://dipterix.org/ravetools/reference/vcg_max_edge_length.md),
[`vcg_subdivision`](https://dipterix.org/ravetools/reference/vcg_subdivision.md)

## Examples

``` r
if (is_not_cran()) {

  sphere <- vcg_sphere()
  cur_max <- vcg_max_edge_length(sphere)
  sphere2 <- vcg_subdivide_max_edge_length(sphere, max_edge_len = cur_max * 0.4)
  vcg_max_edge_length(sphere2)  # should be <= cur_max * 0.4

}
#> [1] 0.04116183
```
