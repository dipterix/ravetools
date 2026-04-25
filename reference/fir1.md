# Window-based `FIR` filter design

Generate a `fir1` filter that is checked against `Matlab` `fir1`
function.

## Usage

``` r
fir1(
  n,
  w,
  type = c("low", "high", "stop", "pass", "DC-0", "DC-1"),
  window = hamming,
  scale = TRUE,
  hilbert = FALSE
)
```

## Arguments

- n:

  filter order

- w:

  band edges, non-decreasing vector in the range 0 to 1, where 1 is the
  `Nyquist` frequency. A scalar for high-pass or low-pass filters, a
  vector pair for band-pass or band-stop, or a vector for an alternating
  pass/stop filter.

- type:

  type of the filter, one of `"low"` for a low-pass filter, `"high"` for
  a high-pass filter, `"stop"` for a stop-band (band-reject) filter,
  `"pass"` for a pass-band filter, `"DC-0"` for a band-pass as the first
  band of a multi-band filter, or `"DC-1"` for a band-stop as the first
  band of a multi-band filter; default `"low"`

- window:

  smoothing window function or a numerical vector. The filter is the
  same shape as the smoothing window. When `window` is a function,
  `window(n+1)` will be called, otherwise the length of the window
  vector needs to have length of `n+1`; default: `hamming`

- scale:

  whether to scale the filter; default is true

- hilbert:

  whether to use 'Hilbert' transformer; default is false

## Value

The `FIR` filter coefficients with class `'Arma'`. The moving average
coefficient is a vector of length `n+1`.
