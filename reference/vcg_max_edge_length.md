# Maximum edge length of a triangular mesh

Returns the length of the longest edge in the mesh.

## Usage

``` r
vcg_max_edge_length(mesh)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`.

## Value

A single numeric value: the maximum edge length, in mesh units.

## See also

[`vcg_average_edge_length`](https://dipterix.org/ravetools/reference/vcg_average_edge_length.md),
[`vcg_subdivide_max_edge_length`](https://dipterix.org/ravetools/reference/vcg_subdivide_max_edge_length.md)

## Examples

``` r
if (is_not_cran()) {

  sphere <- vcg_sphere()
  vcg_max_edge_length(sphere)

}
#> [1] 0.1646472
```
