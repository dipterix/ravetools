# Count boundary and non-manifold edges of a triangular mesh

Detects topology defects that prevent a mesh from being a closed,
manifold, genus-0 surface, a hard precondition of algorithms such as
[`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md).
An edge is a *boundary* edge when it is referenced by exactly one face
(i.e. it bounds a hole), and *non-manifold* when it is referenced by
more than two faces.

## Usage

``` r
vcg_count_edge_defects(mesh)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`.

## Value

A named list with elements `boundary_edges` (number of boundary edges),
`nonmanifold_edges` (number of non-manifold edges), and
`is_closed_manifold` (`TRUE` when both counts are zero, i.e. the mesh is
closed and manifold and ready for
[`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md)).

## Examples

``` r
if (is_not_cran()) {

  sphere <- vcg_sphere()
  vcg_count_edge_defects(sphere)

  defective <- vcg_isosurface(left_hippocampus_mask)
  vcg_count_edge_defects(defective)

}
#> $boundary_edges
#> [1] 40
#> 
#> $nonmanifold_edges
#> [1] 0
#> 
#> $is_closed_manifold
#> [1] FALSE
#> 
```
