# Create a two-dimensional plane in three dimensional space

Create a two-dimensional plane in three dimensional space

## Usage

``` r
plane_geometry(width = 1, height = 1, shape = c(2, 2))
```

## Arguments

- width, height:

  width and height of the plane, must not be `NA`

- shape:

  length of two to indicate the number of vertices along width and
  height, default is only `c(2, 2)` (2 vertices each side, hence one
  grid)

## Value

A triangular mesh of class `'mesh3d'`

## Examples

``` r

plane <- plane_geometry(5, 10, c(12, 22))

if(FALSE) {

  rgl_view({

    rgl_call("shade3d", plane, col = 3)
    rgl_call("wire3d", plane, col = 1)

  })

}
```
