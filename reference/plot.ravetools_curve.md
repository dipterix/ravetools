# Plot method for `ravetools_curve`

Plots a `ravetools_curve` object created by
[`catmull_rom_3d`](https://dipterix.org/ravetools/reference/catmull_rom_3d.md).
When the `rgl` package is available and `use_rgl = TRUE` (default), an
interactive 3D scene is opened. Otherwise three 2-D projection panels
(*x-y*, *x-z*, *y-z*) are drawn using base R graphics.

## Usage

``` r
# S3 method for class 'ravetools_curve'
plot(x, n = 200L, col = "steelblue", pch = 19L, cex = 1, use_rgl = TRUE, ...)
```

## Arguments

- x:

  an object of class `ravetools_curve`.

- n:

  integer; number of sample points used to draw the smooth curve.
  Default is `200L`.

- col:

  color for the spline curve line. Default `"steelblue"`.

- pch, cex:

  plotting character and scaling for the key control points (base-R
  fallback only). Default `pch = 19`, `cex = 1`.

- use_rgl:

  logical; if `TRUE` (default) and the `rgl` package is installed, an
  interactive 3D window is used. Set to `FALSE` to force the base-R 2D
  projection panels.

- ...:

  additional graphical parameters forwarded to the underlying plot
  calls.

## Value

Invisibly returns `x`.
