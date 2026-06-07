# Map continuous values to colors

Linearly maps a numeric vector onto a color ramp, clamping values
outside the given range to the range's endpoints.

## Usage

``` r
color_ramp_continuous(
  values,
  clim = range(values, na.rm = TRUE),
  cmap = grDevices::hcl.colors(11),
  ...
)
```

## Arguments

- values:

  numeric vector of values to map to colors

- clim:

  length-two numeric vector giving the value range to map from; values
  outside `[clim[1], clim[2]]` are clamped to the nearer endpoint before
  mapping. Default is `range(values, na.rm = TRUE)`

- cmap:

  the color ramp to map onto: either a vector of colors (passed to
  [`colorRamp`](https://rdrr.io/r/grDevices/colorRamp.html) to build the
  ramp function) or a function such as the one returned by
  [`colorRamp`](https://rdrr.io/r/grDevices/colorRamp.html) that takes a
  numeric vector with elements in `[0, 1]` and returns an `n x 3` (or
  `n x 4`, with alpha) matrix of color-channel values in `[0, 255]`.
  Default is `grDevices::hcl.colors(11)`

- ...:

  passed to [`colorRamp`](https://rdrr.io/r/grDevices/colorRamp.html)
  when `cmap` is a vector of colors (for example, `alpha = TRUE` to
  build an alpha-aware ramp)

## Value

A character vector of `'#RRGGBB'` (or `'#RRGGBBAA'`) color strings, the
same length as `values`.
