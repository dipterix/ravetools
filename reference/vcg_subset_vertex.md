# Subset mesh by vertex

Subset mesh by vertex

## Usage

``` r
vcg_subset_vertex(x, selector)
```

## Arguments

- x:

  surface mesh

- selector:

  logical vector (must not contain NA), and length must be consistent
  with the number of vertices in `x`: which nodes are to be kept

## Value

A triangular mesh of class `'mesh3d'`, a subset of `x`

## Examples

``` r

sphere <- vcg_sphere()

nv <- ncol(sphere$vb)

selector <- seq_len(nv) > (nv / 2)

sub <- vcg_subset_vertex(sphere, selector)

if(is_not_cran()) {
  rgl_view({

    # subset sphere will be displayed in red
    rgl_call("shade3d", sub, col = 'red')

    # Original sphere will be displayed as wireframe
    rgl_call("wire3d", sphere, col = (2 - selector))

  })
}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’

```
