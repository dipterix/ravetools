# Correct the intensity non-uniformity (bias field) of a 3D volume with the `N4` algorithm

Native re-implementation of the `N4` bias-field correction (`Tustison`
and co-authors, 2010), wrapped the way the `abp_n4` function of 'ANTsPy'
is: outlier intensities are truncated at histogram quantiles, then the
smooth multiplicative bias field is estimated in the log domain by
alternating histogram sharpening with a multi-resolution `B-spline` fit.
The defaults reproduce `abp_n4` of 'ANTsPy' `0.6.3` (the whole image
drives the fit, with one spline span per axis); `mask = "auto"` together
with `spline_distance = 200` gives the older recipe of 'ANTsR' `abpN4`
and 'ANTsPy' `0.3` (an automatic head mask and a 200-millimeter spline
distance). No external dependency is needed.

## Usage

``` r
bias_correction_n4(
  volume,
  vox2ras = NULL,
  mask = NULL,
  weight_mask = NULL,
  intensity_truncation = c(0.025, 0.975, 256),
  shrink_factor = 4,
  iterations = c(50, 50, 50, 50),
  tolerance = 1e-07,
  spline_distance = NULL,
  spline_order = 3,
  histogram_bins = 200,
  bias_fwhm = 0.15,
  wiener_noise = 0.01,
  rescale_intensities = FALSE,
  return_bias_field = FALSE,
  verbose = FALSE
)
```

## Arguments

- volume:

  a 3D numeric array (for example a `'T1'`-weighted `'MRI'`); integer
  and logical arrays are converted to double

- vox2ras:

  optional `4x4` (or `3x4`) matrix mapping the 0-indexed voxel index to
  the anatomical `RAS` coordinate system; if `NULL`, the `"vox2ras"`
  attribute of `volume` is used when present. Only the voxel spacing
  (the norms of its three columns) matters, so that `spline_distance` is
  honored in millimeters; without any geometry the spacing is taken to
  be 1 unit

- mask:

  which voxels drive the bias-field estimate: `NULL` (default) uses the
  whole image, as 'ANTsPy' `0.6.3` does; `"auto"` derives a head mask
  from the (truncated) image the way 'ANTsPy' `get_mask` does (threshold
  at the image mean, erode by two voxels, keep the largest connected
  component, dilate by two voxels, and fill holes, with less clean-up if
  that leaves nothing), which was the default of 'ANTsPy' `0.3`; or an
  array of the same dimensions as `volume` whose non-zero (`TRUE`)
  voxels drive the fit

- weight_mask:

  optional non-negative array of the same dimensions as `volume` giving
  a per-voxel confidence for the `B-spline` fit (for instance a
  white-matter probability map); voxels with zero weight are ignored,
  exactly as if they were outside `mask`

- intensity_truncation:

  numeric vector `c(lower_quantile, upper_quantile, bins)` (default
  `c(0.025, 0.975, 256)`): before anything else, the intensities are
  clamped to the two quantiles of a `bins`-bin histogram of the finite
  voxels, as the `TruncateIntensity` operation of 'ANTsPy' `0.6.3` does
  (see 'Details'); use `NULL` to skip the truncation

- shrink_factor:

  integer (default `4`); the bias field is estimated on the image
  sub-sampled by this factor along every axis (the bias is smooth, so
  this mostly saves time), then evaluated at full resolution

- iterations:

  integer vector; one entry per fitting level giving the maximum number
  of iterations at that level (default `c(50, 50, 50, 50)`, i.e. four
  levels). The `B-spline` control-point mesh doubles its resolution from
  one level to the next

- tolerance:

  convergence threshold (default `1e-7`): a level stops early once the
  coefficient of variation of the ratio between two successive
  bias-field estimates (inside the mask) drops to this value or below

- spline_distance:

  distance between `B-spline` control points at the first level, in the
  units of `vox2ras` (millimeters); a single number or one value per
  axis. The default `NULL` places exactly one spline span over the image
  itself along each axis, without padding (the `spline_param = NULL`
  default of 'ANTsPy' `0.6.3`, that is `-b [1x1x1]`, four cubic control
  points per axis at the first level), which is the same as
  `(dim(volume) - 1) * spacing`. A number such as `200` (the `-b [200]`
  setting of the 'ANTs' `N4BiasFieldCorrection` program and the 'ANTsPy'
  `0.3` default) instead pads the image virtually so that its extent is
  a whole number of spline spans (the 'ANTs' padding rule; a distance
  that divides the field of view exactly needs no padding). The padded
  grid is limited to `2^31 - 1` voxels per axis and the control-point
  lattice at the finest level to `2^27` points in total; values beyond
  these limits (a distance far larger than the field of view, or far
  smaller than the voxel size) raise an error, and a lattice with more
  than eight control points per sample of the sub-sampled fitting grid
  raises a warning, since it usually means the distance was given in the
  wrong units

- spline_order:

  degree of the `B-spline`: `1`, `2` or `3` (cubic, default)

- histogram_bins:

  number of bins of the log-intensity histogram used by the sharpening
  step (default `200`)

- bias_fwhm:

  full width at half maximum, in log-intensity units, of the Gaussian
  that models the bias-field blurring of the histogram (default `0.15`)

- wiener_noise:

  noise constant of the `Wiener` `deconvolution` filter used to sharpen
  the histogram (default `0.01`)

- rescale_intensities:

  logical (default `FALSE`); if `TRUE`, the corrected intensities inside
  the mask are linearly mapped back onto the intensity range that the
  (truncated) input had inside the mask

- return_bias_field:

  logical (default `FALSE`); if `TRUE`, the estimated multiplicative
  bias field is returned **instead of** the corrected volume (the
  'ANTsPy' convention)

- verbose:

  logical (default `FALSE`); print the grid geometry and the convergence
  value of every iteration

## Value

A 3D double array with the same dimensions as `volume`: the
bias-corrected (and truncated) image, or the multiplicative bias field
when `return_bias_field = TRUE`. The `"vox2ras"` attribute is set
whenever a geometry is known (from the argument or from the input
attribute).

## Details

The processing steps are, in order:

1.  **Truncation.** Unless `intensity_truncation` is `NULL`, a histogram
    of the finite voxels is built with `bins` equal-width bins spanning
    their minimum and maximum, except that a minimum of exactly zero is
    raised to `1e-6`, so that an exact-zero background (and anything
    below `1e-6`) is left out; the two requested quantiles are read off
    this histogram by linear interpolation inside the bin, and the
    *whole* image is clamped to that range (no rescaling). This is the
    rule of the installed 'ANTsPy' `0.6.3`; for a non-negative image
    whose background is exactly zero it lifts the background to the
    lower bound, so that every voxel is positive afterwards.

2.  **Mask.** By default (`mask = NULL`) every voxel is in the mask.
    With `mask = "auto"` the mask is derived from the truncated image:
    voxels at or above the image mean, eroded with a ball of radius two,
    reduced to the largest face-connected component, dilated back with
    the same ball and hole-filled. If that yields an empty (or full)
    mask, the clean-up is retried with radius one and then without any
    clean-up.

3.  **`N4`.** The image is sub-sampled by `shrink_factor` and
    transformed to the log domain. At every iteration the histogram of
    the current estimate of the bias-free log image (inside the mask,
    with `histogram_bins` bins) is sharpened by `Wiener` `deconvolution`
    with a Gaussian of width `bias_fwhm`, the expected bias-free
    intensity of each voxel is computed from the sharpened histogram,
    and the difference between the current image and that expectation
    (the residual log bias) is approximated by a `B-spline` field (the
    scattered-data approximation of Lee, `Wolberg` and Shin, 1997,
    weighted by `weight_mask`) that is added to the running
    control-point lattice. A level ends after `iterations[level]`
    iterations or when the convergence measure drops to `tolerance`; the
    lattice is then refined (its number of spans doubles) for the next
    level. The `B-spline` domain is the image itself (default, one span
    per axis) or, for a numeric `spline_distance`, the image padded so
    that its extent is a whole number of spans, which matches the 'ANTs'
    `N4BiasFieldCorrection` program; both the fit on the sub-sampled
    grid and the final evaluation on the full grid use that same domain.

4.  **Output.** The log bias field of the final lattice is evaluated at
    every voxel of the original grid and exponentiated; the corrected
    image is the (truncated) input divided by it everywhere, including
    outside the mask ('ANTs' leaves voxels outside the mask untouched,
    which only differs for background voxels).

Every finite voxel inside the mask (with a positive weight) drives the
fit. As in the 'ITK' filter behind 'ANTs', a positive voxel enters
through its logarithm while a zero or negative one enters with its raw
value; such voxels remain after the default truncation only when the
image has negative intensities (or when `intensity_truncation = NULL`),
and an explicit `mask` (or `mask = "auto"`) keeps them out. Voxels that
are not finite never drive the fit. Every voxel, inside or outside the
mask, is divided by the bias field. Results are deterministic and
identical for any number of threads
([`ravetools_threads`](https://dipterix.org/ravetools/reference/parallel-options.md)).

**Relation to 'ANTs'.** The default arguments (truncation at the `0.025`
and `0.975` quantiles of a 256-bin histogram, the whole image as mask,
shrink factor 4, four levels of at most 50 iterations, tolerance `1e-7`
and one spline span per axis over the image itself) are those of
`abp_n4` in 'ANTsPy' `0.6.3`, which runs `n4_bias_field_correction` with
`mask = NULL` (the whole image) and `spline_param = NULL`
(`-b [1x1x1]`). The `abpN4` function of 'ANTsR' and `abp_n4` of 'ANTsPy'
up to version `0.3` used a `get_mask` head mask and a 200-millimeter
spline distance instead; `mask = "auto", spline_distance = 200`
reproduces that recipe (its truncation then still follows the `0.6.3`
rule above, which differs from the older rule only for images with
negative intensities or values between zero and `1e-6`). The remaining
differences from 'ANTs' are numerical: 'ANTs' computes in single
precision, and with a shrink factor above one it fits on the sub-sampled
image's own, slightly smaller domain, while this function uses the full
image domain at every level.

## References

`Tustison`, N. J., `Avants`, B. B., Cook, P. A., `Zheng`, Y., `Egan`,
A., `Yushkevich`, P. A. and Gee, J. C. (2010). `N4ITK`: improved `N3`
bias correction. *IEEE Transactions on Medical Imaging*, 29(6),
1310-1320.
[doi:10.1109/TMI.2010.2046908](https://doi.org/10.1109/TMI.2010.2046908)

Sled, J. G., `Zijdenbos`, A. P. and Evans, A. C. (1998). A nonparametric
method for automatic correction of intensity `nonuniformity` in `'MRI'`
data. *IEEE Transactions on Medical Imaging*, 17(1), 87-97.
[doi:10.1109/42.668698](https://doi.org/10.1109/42.668698)

Lee, S., `Wolberg`, G. and Shin, S. Y. (1997). Scattered data
interpolation with multilevel `B-splines`. *IEEE Transactions on
Visualization and Computer Graphics*, 3(3), 228-244.
[doi:10.1109/2945.817351](https://doi.org/10.1109/2945.817351)

## See also

[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md)

## Examples

``` r

# Toy phantom: two tissue classes inside a sphere, multiplied by a
# smooth bias that increases along x
nd <- c(24, 24, 24)
g <- expand.grid(x = 0:23, y = 0:23, z = 0:23)
r <- sqrt((g$x - 11.5)^2 + (g$y - 11.5)^2 + (g$z - 11.5)^2)
tissue <- ifelse(r < 6, 150, ifelse(r < 10, 80, 0))
bias <- exp(0.4 * (g$x - 11.5) / 24)
set.seed(1)
volume <- array(tissue * bias + rnorm(nrow(g), sd = 2), nd)
vox2ras <- diag(c(2, 2, 2, 1))        # 2 mm isotropic voxels

# defaults of 'ANTsPy' 0.6.3 abp_n4: the whole image drives the fit,
# with one spline span per axis (fewer iterations to keep this fast)
corrected <- bias_correction_n4(
  volume, vox2ras = vox2ras, shrink_factor = 2, iterations = c(20, 20))

# the estimated field follows the true bias inside the head
est <- bias_correction_n4(
  volume, vox2ras = vox2ras, return_bias_field = TRUE,
  shrink_factor = 2, iterations = c(20, 20))
inside <- r < 10
cor(log(est[inside]), log(bias[inside]))
#> [1] 0.9918538

# the corrected tissue intensities are flatter than the input
sd(volume[r < 6]) > sd(corrected[r < 6])
#> [1] TRUE

# the older 'ANTsR' / 'ANTsPy' 0.3 recipe: an automatic head mask and a
# spline distance in millimeters (40 mm suits this small toy; 200 mm is
# the usual value for a head)
est2 <- bias_correction_n4(
  volume, vox2ras = vox2ras, mask = "auto", spline_distance = 40,
  return_bias_field = TRUE, shrink_factor = 2, iterations = c(20, 20))
cor(log(est2[inside]), log(bias[inside]))
#> [1] 0.9985041
```
