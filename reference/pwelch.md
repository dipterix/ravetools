# Calculate 'Welch Periodogram'

`pwelch` is for single signal trace only; `mv_pwelch` is for multiple
traces. Currently `mv_pwelch` is experimental and should not be called
directly.

## Usage

``` r
pwelch(
  x,
  fs,
  window = 64,
  noverlap = window/2,
  nfft = "auto",
  window_family = hamming,
  col = "black",
  xlim = NULL,
  ylim = NULL,
  main = "Welch periodogram",
  plot = 0,
  log = c("xy", "", "x", "y"),
  ...
)

# S3 method for class '`ravetools-pwelch`'
print(x, ...)

# S3 method for class '`ravetools-pwelch`'
plot(
  x,
  log = c("xy", "x", "y", ""),
  se = FALSE,
  xticks,
  type = "l",
  add = FALSE,
  col = graphics::par("fg"),
  col.se = "orange",
  alpha.se = 0.5,
  lty = 1,
  lwd = 1,
  cex = 1,
  las = 1,
  main = "Welch periodogram",
  xlab,
  ylab,
  xlim = NULL,
  ylim = NULL,
  xaxs = "i",
  yaxs = "i",
  xline = 1.2 * cex,
  yline = 2 * cex,
  mar = c(2.6, 3.8, 2.1, 0.6) * (0.5 + cex/2),
  mgp = cex * c(2, 0.5, 0),
  tck = -0.02 * cex,
  grid = TRUE,
  ...
)

mv_pwelch(
  x,
  margin,
  fs,
  window = 64,
  noverlap = window/2,
  nfft = "auto",
  window_family = hamming
)
```

## Arguments

- x:

  numerical vector or a row-major vector, signals. If `x` is a matrix,
  then each row is a channel. For `plot` function, `x` is the instance
  returned by `pwelch` function.

- fs:

  sample rate, average number of time points per second

- window:

  window length in time points, default size is `64`

- noverlap:

  overlap between two adjacent windows, measured in time points; default
  is half of the `window`

- nfft:

  number of points in window function; default is automatically
  determined from input data and window, scaled up to the nearest power
  of 2

- window_family:

  function generator for generating filter windows, default is
  [`hamming`](https://dipterix.org/ravetools/reference/filter-window.md).
  This can be any window function listed in the filter window family, or
  any window generator function from package `gsignal`. Default is
  [`hamming`](https://dipterix.org/ravetools/reference/filter-window.md).
  For 'iEEG' users, both `hamming` and
  [`blackmanharris`](https://dipterix.org/ravetools/reference/filter-window.md)
  are offered by 'EEG-lab'; while `blackmanharris` offers better
  attenuation than Hamming windows, it also has lower spectral
  resolution. `hamming` has a 42.5 dB side-lobe attenuation. This may
  mask spectral content below this value (relative to the peak spectral
  content). Choosing different windows enables you to make trade-off
  between resolution (e.g., using a rectangular window) and side-lobe
  attenuation (e.g., using a
  [`hanning`](https://dipterix.org/ravetools/reference/filter-window.md)
  window)

- col, xlim, ylim, main, type, cex, las, xlab, ylab, lty, lwd, xaxs,
  yaxs, mar, mgp, tck:

  parameters passed to
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html)

- plot:

  integer, whether to plot the result or not; choices are `0`, no plot;
  `1` plot on a new canvas; `2` add to existing canvas

- log:

  indicates which axis should be `log10`-transformed, used by the plot
  function. For `'x'` axis, it's `log10`-transform; for `'y'` axis, it's
  `10log10`-transform (decibel unit). Choices are `"xy"`, `"x"`, `"y"`,
  and `""`.

- ...:

  will be passed to `plot.pwelch` or ignored

- se:

  logical or a positive number indicating whether to plot standard error
  of mean; default is false. If provided with a number, then a multiple
  of standard error will be drawn. This option is only available when
  power is in log-scale (decibel unit)

- xticks:

  ticks to show on frequency axis

- add:

  logical, whether the plot should be added to existing canvas

- col.se, alpha.se:

  controls the color and opacity of the standard error

- xline, yline:

  controls how close the axis labels to the corresponding axes

- grid:

  whether to draw rectangular grid lines to the plot; only respected
  when `add=FALSE`; default is true

- margin:

  the margin in which `pwelch` should be applied to

## Value

A list with class `'ravetools-pwelch'` that contains the following
items:

- `freq`:

  frequencies used to calculate the 'periodogram'

- `spec`:

  resulting spectral power for each frequency

- `window`:

  window function(in numerical vector) used

- `noverlap`:

  number of overlapping time-points between two adjacent windows

- `nfft`:

  number of basis functions

- `fs`:

  sample rate

- `x_len`:

  input signal length

- `method`:

  a character string `'Welch'`

## Examples

``` r
x <- rnorm(1000)
pwel <- pwelch(x, 100)
pwel
#> Welch Periodogram:
#>   # channels: 1
#>   time points: 1000
#>   sample rate: 100.00
#>   window size: 64
#>   window overlaps: 32
#>   filter count: 256

plot(pwel, log = "xy")

```
