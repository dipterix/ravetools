# Detect and repair defects in a triangular surface mesh

Repairs common defects that prevent a mesh from being a closed,
manifold, genus-0 surface - a hard precondition of algorithms such as
[`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md).
Typical sources of such defects are surfaces extracted from volumes via
marching-cubes-style algorithms (e.g.
[`vcg_isosurface`](https://dipterix.org/ravetools/reference/vcg_isosurface.md)),
which can leave behind small "cracks": isolated boundary-edge loops
bounding tiny holes that are not closed by simple vertex-welding.

## Usage

``` r
vcg_fix_defects(
  mesh,
  merge_tolerance = NA,
  max_hole_size = 100L,
  verbose = FALSE
)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`.

- merge_tolerance:

  distance (in mesh units) below which vertices are welded together;
  default is `NA`, in which case the tolerance is derived automatically
  as `1e-4` times the mesh's average edge length.

- max_hole_size:

  maximum number of boundary edges of a hole that will be triangulated
  (ear-cutting fill); holes larger than this threshold are left
  untouched (and will be reported as remaining boundary edges in
  `attr(..., "info")`). Default is `100`.

- verbose:

  whether to print a short before/after diagnostic report; default is
  `FALSE`.

## Value

A repaired triangular mesh of class `'mesh3d'`, with an additional
attribute `"info"`, a named list reporting what was found and changed:
`boundary_edges_before/after`, `nonmanifold_edges_before/after`,
`vertices_merged`, `merge_tolerance`, `holes_filled`, `is_oriented`,
`is_orientable`, `normals_flipped_outward`, and `is_closed_manifold`
(`TRUE` when the repaired mesh is closed and manifold, i.e. ready for
[`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md)).

## Details

The repair pipeline applies, in order:

1.  Remove degenerate and duplicate faces.

2.  Weld near-coincident vertices (closes cracks caused by duplicated
    vertices), using `merge_tolerance` or, by default, a distance
    derived from the mesh's average edge length.

3.  Triangulate ("ear-cut fill") any remaining small boundary loops –
    i.e. *isolated edges* / genuine small holes that welding alone
    cannot close, up to `max_hole_size` edges.

4.  Remove unreferenced vertices.

5.  Re-orient all faces coherently (consistent winding order), and, if
    the result is a single watertight component, flip normals to point
    outward (this last step assumes the geometry is meant to be
    watertight).

## Examples

``` r
if (is_not_cran()) {

  mesh <- vcg_isosurface(left_hippocampus_mask)

  repaired <- vcg_fix_defects(mesh, verbose = TRUE)

  attr(repaired, "info")$is_closed_manifold

  # repaired mesh can now be inflated
  inflated <- mris_inflate(repaired, scale_brain = FALSE)

}
#> vcgFixDefects: input nv=4318 nf=8652, boundary edges=40, non-manifold edges=0
#> vcgFixDefects: [1] removed degenerate/duplicate faces -> nv=4318 nf=8652
#> vcgFixDefects: [2] merged 0 close vertices -> nv=4318 nf=8652
#> vcgFixDefects: [3a] topology/normals ready, starting hole fill (max_hole_size=100)
#> vcgFixDefects: [3b] filled 10 hole(s) -> nv=4318 nf=8672
#> vcgFixDefects: [4] removed unreferenced vertices -> nv=4318 nf=8672
#> vcgFixDefects: [5a] topology rebuilt, orienting coherently
#> vcgFixDefects: [5b] oriented=1 orientable=1
#> vcgFixDefects: removed/merged 0 vertices (tol=9.97168e-05), filled 10 hole(s)
#> vcgFixDefects: output nv=4318 nf=8672, boundary edges=0, non-manifold edges=0, oriented=yes, orientable=yes, normals_flipped_outward=no
```
