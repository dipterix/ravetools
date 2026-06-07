# Compute the average edge length of a triangular mesh

Computes the average length of all face edges (each edge is counted once
per incident face, so edges shared by two faces are counted twice).
Useful as a scale-aware reference length, e.g. to derive a vertex-merge
tolerance such as the one used internally by
[`vcg_fix_defects`](https://dipterix.org/ravetools/reference/vcg_fix_defects.md).

## Usage

``` r
vcg_average_edge_length(mesh)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`.

## Value

A single numeric value: the average edge length, in mesh units.

## Examples

``` r
if (is_not_cran()) {

  sphere <- vcg_sphere()
  vcg_average_edge_length(sphere)

}
#> [1] 0.1507297
```
