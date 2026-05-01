# Design an 'IIR' filter

Design an 'IIR' filter

## Usage

``` r
design_filter_iir(
  method = c("butter", "cheby1", "cheby2", "ellip"),
  sample_rate,
  filter_order = NA,
  high_pass_freq = NA,
  high_pass_trans_freq = NA,
  low_pass_freq = NA,
  low_pass_trans_freq = NA,
  passband_ripple = 0.1,
  stopband_attenuation = 40
)
```

## Arguments

- method:

  filter method name, choices are `"butter"`, `"cheby1"`, `"cheby2"`,
  and `"ellip"`

- sample_rate:

  sampling frequency

- filter_order:

  suggested filter order. Notice filters with higher orders may become
  numerically unstable, hence this number is only a suggested number. If
  the filter is unstable, this function will choose a lower order; leave
  this input `NA` (default) if undecided.

- high_pass_freq:

  high-pass frequency; default is `NA` (no high-pass filter will be
  applied)

- high_pass_trans_freq:

  high-pass frequency band-width; default is automatically inferred from
  filter type.

- low_pass_freq:

  low-pass frequency; default is `NA` (no low-pass filter will be
  applied)

- low_pass_trans_freq:

  low-pass frequency band-width; default is automatically inferred from
  filter type.

- passband_ripple:

  allowable pass-band ripple in decibel; default is `0.1`

- stopband_attenuation:

  minimum stop-band attenuation (in decibel) at transition frequency;
  default is `40` dB.

## Value

A filter in 'Arma' form.

## Examples

``` r

sample_rate <- 500

my_diagnose <- function(
    filter, vlines = c(8, 12), cutoffs = c(-3, -6)) {
  diagnose_filter(
    b = filter$b,
    a = filter$a,
    fs = sample_rate,
    vlines = vlines,
    cutoffs = cutoffs
  )
}

# ---- Default using butterworth to generate 8-12 bandpass filter ----

# Butterworth filter with cut-off frequency
# 7 ~ 13 (default transition bandwidth is 1Hz) at -3 dB
filter <- design_filter_iir(
  method = "butter",
  low_pass_freq = 12,
  high_pass_freq = 8,
  sample_rate = 500
)

filter
#> <IIR filter>
#>   Type : pass
#>   Method: butter
#>   Order: 5
#>   Magnitude:
#>     Freq=8 Hz, mag=-0.1059 dB (expected=-0.1 dB)
#>     Freq=12 Hz, mag=-0.1019 dB (expected=-0.1 dB)
#>     Freq=6 Hz, mag=-23.46 dB (expected=-40 dB)
#>     Freq=15 Hz, mag=-17.06 dB (expected=-40 dB)
#>   Reciprocal condition number: 1.8e-15 > .Machine$double.eps

my_diagnose(filter)


## explicit bandwidths and attenuation (sharper transition)

# Butterworth filter with cut-off frequency
# passband ripple is 0.5 dB (8-12 Hz)
# stopband attenuation is 40 dB (5-18 Hz)
filter <- design_filter_iir(
  method = "butter",
  low_pass_freq = 12, low_pass_trans_freq = 6,
  high_pass_freq = 8, high_pass_trans_freq = 3,
  sample_rate = 500,
  passband_ripple = 0.5,
  stopband_attenuation = 40
)

filter
#> <IIR filter>
#>   Type : pass
#>   Method: butter
#>   Order: 5
#>   Magnitude:
#>     Freq=8 Hz, mag=-0.5076 dB (expected=-0.5 dB)
#>     Freq=12 Hz, mag=-0.5011 dB (expected=-0.5 dB)
#>     Freq=5 Hz, mag=-45.85 dB (expected=-40 dB)
#>     Freq=18 Hz, mag=-41.04 dB (expected=-40 dB)
#>   Reciprocal condition number: 1.7e-15 > .Machine$double.eps

my_diagnose(filter)


# ---- cheby1 --------------------------------

filter <- design_filter_iir(
  method = "cheby1",
  low_pass_freq = 12,
  high_pass_freq = 8,
  sample_rate = 500
)

my_diagnose(filter)


# ---- cheby2 --------------------------------

filter <- design_filter_iir(
  method = "cheby2",
  low_pass_freq = 12,
  high_pass_freq = 8,
  sample_rate = 500
)

my_diagnose(filter)


# ----- ellip ---------------------------------

filter <- design_filter_iir(
  method = "ellip",
  low_pass_freq = 12,
  high_pass_freq = 8,
  sample_rate = 500
)

my_diagnose(filter)




```
