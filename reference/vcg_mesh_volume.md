# Compute volume for manifold meshes

Compute volume for manifold meshes

## Usage

``` r
vcg_mesh_volume(mesh)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`

## Value

The numeric volume of the mesh

## Examples

``` r
# Initial mesh
mesh <- vcg_sphere()

vcg_mesh_volume(mesh)
#> [1] 4.152741
```
