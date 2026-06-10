# Simple 3-dimensional sphere mesh

Simple 3-dimensional sphere mesh

## Usage

``` r
vcg_sphere(sub_division = 3L, normals = TRUE)
```

## Arguments

- sub_division:

  density of vertex in the resulting mesh

- normals:

  whether the normal vectors should be calculated

## Value

A `'mesh3d'` object

## Examples

``` r

vcg_sphere()
#>  mesh3d object with 642 vertices, 1280 triangles.
```
