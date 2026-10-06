# N4 bias-field correction with the 'ANTsPy' `abp_n4` preprocessing (intensity
# truncation, then N4). The defaults follow the installed ANTsPy 0.6.3
# (whole-image mask, one B-spline span per axis); mask = "auto" plus
# spline_distance = 200 gives the older ANTsR / ANTsPy 0.3 recipe. Native C++
# core in src/imaging_intensity_n4.cpp; this file validates arguments, builds
# the mask, and assembles the corrected volume / bias field.

#' @title Correct the intensity non-uniformity (bias field) of a 3D volume
#' with the \code{N4} algorithm
#' @description
#' Native re-implementation of the \code{N4} bias-field correction
#' (\verb{Tustison} and co-authors, 2010), wrapped the way the \code{abp_n4}
#' function of \pkg{'ANTsPy'} is: outlier intensities are truncated at
#' histogram quantiles, then the smooth multiplicative bias field is estimated
#' in the log domain by alternating histogram sharpening with a
#' multi-resolution \verb{B-spline} fit. The defaults reproduce
#' \code{abp_n4} of \pkg{'ANTsPy'} \code{0.6.3} (the whole image drives the fit,
#' with one spline span per axis); \code{mask = "auto"} together with
#' \code{spline_distance = 200} gives the older recipe of \pkg{'ANTsR'}
#' \code{abpN4} and \pkg{'ANTsPy'} \code{0.3} (an automatic head mask and a
#' 200-millimeter spline distance). No external dependency is needed.
#'
#' @param volume a 3D numeric array (for example a \code{'T1'}-weighted
#' \code{'MRI'}); integer and logical arrays are converted to double
#' @param vox2ras optional \code{4x4} (or \code{3x4}) matrix mapping the
#' 0-indexed voxel index to the anatomical \verb{RAS} coordinate system; if
#' \code{NULL}, the \code{"vox2ras"} attribute of \code{volume} is used when
#' present. Only the voxel spacing (the norms of its three columns) matters, so
#' that \code{spline_distance} is honored in millimeters; without any geometry
#' the spacing is taken to be 1 unit
#' @param mask which voxels drive the bias-field estimate: \code{NULL}
#' (default) uses the whole image, as \pkg{'ANTsPy'} \code{0.6.3} does;
#' \code{"auto"} derives a head mask from the (truncated) image the way
#' \pkg{'ANTsPy'} \code{get_mask} does (threshold at the image mean, erode by
#' two voxels, keep the largest connected component, dilate by two voxels, and
#' fill holes, with less clean-up if that leaves nothing), which was the
#' default of \pkg{'ANTsPy'} \code{0.3}; or an array of the same dimensions as
#' \code{volume} whose non-zero (\code{TRUE}) voxels drive the fit
#' @param weight_mask optional non-negative array of the same dimensions as
#' \code{volume} giving a per-voxel confidence for the \verb{B-spline} fit (for
#' instance a white-matter probability map); voxels with zero weight are
#' ignored, exactly as if they were outside \code{mask}
#' @param intensity_truncation numeric vector
#' \code{c(lower_quantile, upper_quantile, bins)} (default
#' \code{c(0.025, 0.975, 256)}): before anything else, the intensities are
#' clamped to the two quantiles of a \code{bins}-bin histogram of the finite
#' voxels, as the \code{TruncateIntensity} operation of \pkg{'ANTsPy'}
#' \code{0.6.3} does (see 'Details'); use \code{NULL} to skip the truncation
#' @param shrink_factor integer (default \code{4}); the bias field is estimated
#' on the image sub-sampled by this factor along every axis (the bias is smooth,
#' so this mostly saves time), then evaluated at full resolution
#' @param iterations integer vector; one entry per fitting level giving the
#' maximum number of iterations at that level (default \code{c(50, 50, 50,
#' 50)}, i.e. four levels). The \verb{B-spline} control-point mesh doubles its
#' resolution from one level to the next
#' @param tolerance convergence threshold (default \code{1e-7}): a level stops
#' early once the coefficient of variation of the ratio between two successive
#' bias-field estimates (inside the mask) drops to this value or below
#' @param spline_distance distance between \verb{B-spline} control points at
#' the first level, in the units of \code{vox2ras} (millimeters); a single
#' number or one value per axis. The default \code{NULL} places exactly one
#' spline span over the image itself along each axis, without padding (the
#' \code{spline_param = NULL} default of \pkg{'ANTsPy'} \code{0.6.3}, that is
#' \code{-b [1x1x1]}, four cubic control points per axis at the first level),
#' which is the same as \code{(dim(volume) - 1) * spacing}. A number such as
#' \code{200} (the \code{-b [200]} setting of the \pkg{'ANTs'}
#' \code{N4BiasFieldCorrection} program and the \pkg{'ANTsPy'} \code{0.3}
#' default) instead pads the image virtually so that its extent is a whole
#' number of spline spans (the \pkg{'ANTs'} padding rule; a distance that
#' divides the field of view exactly needs no padding). The padded grid is
#' limited to \code{2^31 - 1} voxels per axis and the control-point lattice at
#' the finest level to \code{2^27} points in total; values beyond these limits
#' (a distance far larger than the field of view, or far smaller than the voxel
#' size) raise an error, and a lattice with more than eight control points per
#' sample of the sub-sampled fitting grid raises a warning, since it usually
#' means the distance was given in the wrong units
#' @param spline_order degree of the \verb{B-spline}: \code{1}, \code{2} or
#' \code{3} (cubic, default)
#' @param histogram_bins number of bins of the log-intensity histogram used by
#' the sharpening step (default \code{200})
#' @param bias_fwhm full width at half maximum, in log-intensity units, of the
#' Gaussian that models the bias-field blurring of the histogram (default
#' \code{0.15})
#' @param wiener_noise noise constant of the \verb{Wiener} \verb{deconvolution}
#' filter used to sharpen the histogram (default \code{0.01})
#' @param rescale_intensities logical (default \code{FALSE}); if \code{TRUE},
#' the corrected intensities inside the mask are linearly mapped back onto the
#' intensity range that the (truncated) input had inside the mask
#' @param return_bias_field logical (default \code{FALSE}); if \code{TRUE}, the
#' estimated multiplicative bias field is returned \strong{instead of} the
#' corrected volume (the \pkg{'ANTsPy'} convention)
#' @param verbose logical (default \code{FALSE}); print the grid geometry and
#' the convergence value of every iteration
#'
#' @details
#' The processing steps are, in order:
#' \enumerate{
#' \item \strong{Truncation.} Unless \code{intensity_truncation} is
#' \code{NULL}, a histogram of the finite voxels is built with \code{bins}
#' equal-width bins spanning their minimum and maximum, except that a minimum
#' of exactly zero is raised to \code{1e-6}, so that an exact-zero background
#' (and anything below \code{1e-6}) is left out; the two requested quantiles
#' are read off this histogram by linear interpolation inside the bin, and the
#' \emph{whole} image is clamped to that range (no rescaling). This is the rule
#' of the installed \pkg{'ANTsPy'} \code{0.6.3}; for a non-negative image whose
#' background is exactly zero it lifts the background to the lower bound, so
#' that every voxel is positive afterwards.
#' \item \strong{Mask.} By default (\code{mask = NULL}) every voxel is in the
#' mask. With \code{mask = "auto"} the mask is derived from the truncated
#' image: voxels at or above the image mean, eroded with a ball of radius two,
#' reduced to the largest face-connected component, dilated back with the same
#' ball and hole-filled. If that yields an empty (or full) mask, the clean-up is
#' retried with radius one and then without any clean-up.
#' \item \strong{\code{N4}.} The image is sub-sampled by \code{shrink_factor}
#' and transformed to the log domain. At every iteration the histogram of the
#' current estimate of the bias-free log image (inside the mask, with
#' \code{histogram_bins} bins) is sharpened by \verb{Wiener} \verb{deconvolution} with
#' a Gaussian of width \code{bias_fwhm}, the expected bias-free intensity of
#' each voxel is computed from the sharpened histogram, and the difference
#' between the current image and that expectation (the residual log bias) is
#' approximated by a \verb{B-spline} field (the scattered-data approximation
#' of Lee, \verb{Wolberg} and Shin, 1997, weighted by \code{weight_mask}) that
#' is added to the running control-point lattice. A level ends after
#' \code{iterations[level]} iterations or when the convergence measure drops to
#' \code{tolerance}; the lattice is then refined (its number of spans doubles)
#' for the next level. The \verb{B-spline} domain is the image itself (default,
#' one span per axis) or, for a numeric \code{spline_distance}, the image
#' padded so that its extent is a whole number of spans, which matches the
#' \pkg{'ANTs'} \code{N4BiasFieldCorrection} program; both the fit on the
#' sub-sampled grid and the final evaluation on the full grid use that same
#' domain.
#' \item \strong{Output.} The log bias field of the final lattice is evaluated
#' at every voxel of the original grid and exponentiated; the corrected image
#' is the (truncated) input divided by it everywhere, including outside the
#' mask (\pkg{'ANTs'} leaves voxels outside the mask untouched, which only
#' differs for background voxels).
#' }
#' Every finite voxel inside the mask (with a positive weight) drives the fit.
#' As in the \pkg{'ITK'} filter behind \pkg{'ANTs'}, a positive voxel enters
#' through its logarithm while a zero or negative one enters with its raw
#' value; such voxels remain after the default truncation only when the image
#' has negative intensities (or when \code{intensity_truncation = NULL}), and
#' an explicit \code{mask} (or \code{mask = "auto"}) keeps them out. Voxels
#' that are not finite never drive the fit. Every voxel, inside or outside the
#' mask, is divided by the bias field. Results are deterministic and identical
#' for any number of threads (\code{\link{ravetools_threads}}).
#'
#' \strong{Relation to \pkg{'ANTs'}.} The default arguments (truncation at the
#' \code{0.025} and \code{0.975} quantiles of a 256-bin histogram, the whole
#' image as mask, shrink factor 4, four levels of at most 50 iterations,
#' tolerance \code{1e-7} and one spline span per axis over the image itself)
#' are those of \code{abp_n4} in \pkg{'ANTsPy'} \code{0.6.3}, which runs
#' \code{n4_bias_field_correction} with \code{mask = NULL} (the whole image)
#' and \code{spline_param = NULL} (\code{-b [1x1x1]}). The \code{abpN4}
#' function of \pkg{'ANTsR'} and \code{abp_n4} of \pkg{'ANTsPy'} up to version
#' \code{0.3} used a \code{get_mask} head mask and a 200-millimeter spline
#' distance instead; \code{mask = "auto", spline_distance = 200} reproduces
#' that recipe (its truncation then still follows the \code{0.6.3} rule above,
#' which differs from the older rule only for images with negative intensities
#' or values between zero and \code{1e-6}). The remaining differences from
#' \pkg{'ANTs'} are numerical: \pkg{'ANTs'} computes in single precision, and
#' with a shrink factor above one it fits on the sub-sampled image's own,
#' slightly smaller domain, while this function uses the full image domain at
#' every level.
#'
#' @returns A 3D double array with the same dimensions as \code{volume}: the
#' bias-corrected (and truncated) image, or the multiplicative bias field when
#' \code{return_bias_field = TRUE}. The \code{"vox2ras"} attribute is set
#' whenever a geometry is known (from the argument or from the input
#' attribute).
#'
#' @references
#' \verb{Tustison}, N. J., \verb{Avants}, B. B., Cook, P. A., \verb{Zheng}, Y.,
#' \verb{Egan}, A., \verb{Yushkevich}, P. A. and Gee, J. C. (2010). \verb{N4ITK}:
#' improved \verb{N3} bias correction. \emph{IEEE Transactions on Medical
#' Imaging}, 29(6), 1310-1320. \doi{10.1109/TMI.2010.2046908}
#'
#' Sled, J. G., \verb{Zijdenbos}, A. P. and Evans, A. C. (1998). A
#' nonparametric method for automatic correction of intensity \verb{nonuniformity} in
#' \code{'MRI'} data. \emph{IEEE Transactions on Medical Imaging}, 17(1),
#' 87-97. \doi{10.1109/42.668698}
#'
#' Lee, S., \verb{Wolberg}, G. and Shin, S. Y. (1997). Scattered data
#' interpolation with multilevel \verb{B-splines}. \emph{IEEE Transactions on
#' Visualization and Computer Graphics}, 3(3), 228-244.
#' \doi{10.1109/2945.817351}
#'
#' @seealso \code{\link{register_volume3d}}
#' @examples
#'
#' # Toy phantom: two tissue classes inside a sphere, multiplied by a
#' # smooth bias that increases along x
#' nd <- c(24, 24, 24)
#' g <- expand.grid(x = 0:23, y = 0:23, z = 0:23)
#' r <- sqrt((g$x - 11.5)^2 + (g$y - 11.5)^2 + (g$z - 11.5)^2)
#' tissue <- ifelse(r < 6, 150, ifelse(r < 10, 80, 0))
#' bias <- exp(0.4 * (g$x - 11.5) / 24)
#' set.seed(1)
#' volume <- array(tissue * bias + rnorm(nrow(g), sd = 2), nd)
#' vox2ras <- diag(c(2, 2, 2, 1))        # 2 mm isotropic voxels
#'
#' # defaults of 'ANTsPy' 0.6.3 abp_n4: the whole image drives the fit,
#' # with one spline span per axis (fewer iterations to keep this fast)
#' corrected <- bias_correction_n4(
#'   volume, vox2ras = vox2ras, shrink_factor = 2, iterations = c(20, 20))
#'
#' # the estimated field follows the true bias inside the head
#' est <- bias_correction_n4(
#'   volume, vox2ras = vox2ras, return_bias_field = TRUE,
#'   shrink_factor = 2, iterations = c(20, 20))
#' inside <- r < 10
#' cor(log(est[inside]), log(bias[inside]))
#'
#' # the corrected tissue intensities are flatter than the input
#' sd(volume[r < 6]) > sd(corrected[r < 6])
#'
#' # the older 'ANTsR' / 'ANTsPy' 0.3 recipe: an automatic head mask and a
#' # spline distance in millimeters (40 mm suits this small toy; 200 mm is
#' # the usual value for a head)
#' est2 <- bias_correction_n4(
#'   volume, vox2ras = vox2ras, mask = "auto", spline_distance = 40,
#'   return_bias_field = TRUE, shrink_factor = 2, iterations = c(20, 20))
#' cor(log(est2[inside]), log(bias[inside]))
#'
#' @export
bias_correction_n4 <- function(
    volume, vox2ras = NULL, mask = NULL, weight_mask = NULL,
    intensity_truncation = c(0.025, 0.975, 256), shrink_factor = 4,
    iterations = c(50, 50, 50, 50), tolerance = 1e-7, spline_distance = NULL,
    spline_order = 3, histogram_bins = 200, bias_fwhm = 0.15, wiener_noise = 0.01,
    rescale_intensities = FALSE, return_bias_field = FALSE, verbose = FALSE) {

  fname <- "`bias_correction_n4`"

  # ---- switches ---------------------------------------------------------------
  verbose <- isTRUE(as.logical(verbose))
  rescale_intensities <- isTRUE(as.logical(rescale_intensities))
  return_bias_field <- isTRUE(as.logical(return_bias_field))

  # ---- volume ---------------------------------------------------------------
  if (!typeof(volume) %in% c("double", "integer", "logical")) {
    stop(fname, ": `volume` must be a numeric (or logical) 3D array.")
  }
  d <- dim(volume)
  if (length(d) < 3 || any(d[1:3] < 2) || (length(d) > 3 && !all(d[-(1:3)] == 1))) {
    stop(fname, ": `volume` must be a 3D array with at least 2 voxels along each axis.")
  }
  d <- as.integer(d[1:3])
  n <- prod(d)
  if (is.null(vox2ras)) {
    vox2ras <- attr(volume, "vox2ras")
  }
  x <- as.double(volume)
  if (length(x) != n) {
    stop(fname, ": `volume` must be a 3D array.")
  }

  # ---- geometry -------------------------------------------------------------
  if (is.null(vox2ras)) {
    spacing <- c(1, 1, 1)
  } else {
    vox2ras <- as.matrix(vox2ras)
    if (nrow(vox2ras) == 3L && ncol(vox2ras) == 4L) {
      vox2ras <- rbind(vox2ras, c(0, 0, 0, 1))
    }
    if (!all(dim(vox2ras) == c(4L, 4L)) || !is.numeric(vox2ras) || !all(is.finite(vox2ras))) {
      stop(fname, ": `vox2ras` must be a finite 4x4 (or 3x4) numeric matrix.")
    }
    spacing <- sqrt(colSums(vox2ras[1:3, 1:3, drop = FALSE]^2))
    if (!all(is.finite(spacing)) || any(spacing <= 0)) {
      stop(fname, ": `vox2ras` has a degenerate (zero-length) axis.")
    }
  }

  # ---- scalar settings ------------------------------------------------------
  shrink_factor <- validate_n4_integer(shrink_factor, "shrink_factor", 1L, fname)
  spline_order <- validate_n4_integer(spline_order, "spline_order", 1L, fname)
  if (spline_order > 3L) {
    stop(fname, ": `spline_order` must be 1, 2 or 3.")
  }
  histogram_bins <- validate_n4_integer(histogram_bins, "histogram_bins", 2L, fname)
  if (!is.numeric(iterations) || !length(iterations) || any(!is.finite(iterations)) ||
      any(iterations < 0) || any(iterations != round(iterations))) {
    stop(fname, ": `iterations` must be a vector of non-negative integers (one per level).")
  }
  iterations <- as.integer(iterations)
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) || tolerance < 0) {
    stop(fname, ": `tolerance` must be a single non-negative number.")
  }
  if (is.null(spline_distance)) {
    # one span over the image itself along each axis, without padding: the
    # `spline_param = NULL` (-b [1x1x1]) default of ANTsPy 0.6.3
    spline_distance <- (d - 1) * spacing
  } else if (!is.numeric(spline_distance) || !length(spline_distance) %in% c(1L, 3L) ||
             any(!is.finite(spline_distance)) || any(spline_distance <= 0)) {
    stop(fname, ": `spline_distance` must be NULL or one or three positive numbers.")
  }
  spline_distance <- rep_len(as.double(spline_distance), 3L)
  if (!is.numeric(bias_fwhm) || length(bias_fwhm) != 1L || !is.finite(bias_fwhm) || bias_fwhm <= 0) {
    stop(fname, ": `bias_fwhm` must be a single positive number.")
  }
  if (!is.numeric(wiener_noise) || length(wiener_noise) != 1L || !is.finite(wiener_noise) || wiener_noise <= 0) {
    stop(fname, ": `wiener_noise` must be a single positive number.")
  }

  # ---- 1. truncation --------------------------------------------------------
  if (!is.null(intensity_truncation)) {
    if (!is.numeric(intensity_truncation) || length(intensity_truncation) != 3L ||
        any(!is.finite(intensity_truncation)) ||
        intensity_truncation[1] < 0 || intensity_truncation[2] > 1 ||
        intensity_truncation[1] > intensity_truncation[2] ||
        intensity_truncation[3] < 1 || intensity_truncation[3] != round(intensity_truncation[3])) {
      stop(fname, ": `intensity_truncation` must be c(lower, upper, bins) with 0 <= lower <= upper <= 1 and a whole number of bins >= 1, or NULL.")
    }
    intensity_truncation <- as.double(intensity_truncation)
    bounds <- n4_truncate_quantiles(x, intensity_truncation[1], intensity_truncation[2],
                                    as.integer(intensity_truncation[3]))
    if (all(is.finite(bounds))) {
      lo_idx <- which(x < bounds[1])
      if (length(lo_idx)) x[lo_idx] <- bounds[1]
      hi_idx <- which(x > bounds[2])
      if (length(hi_idx)) x[hi_idx] <- bounds[2]
      if (verbose) {
        message(sprintf("[N4] intensities truncated to [%.6g, %.6g]", bounds[1], bounds[2]))
      }
    }
  }

  # ---- 2. mask and weights ---------------------------------------------------
  if (is.null(mask)) {
    # the whole image drives the fit (ANTsPy 0.6.3 default)
    mask <- rep(TRUE, n)
  } else if (is.character(mask)) {
    if (!identical(mask, "auto")) {
      stop(fname, ": `mask` must be NULL (the whole image), \"auto\" (an automatic ",
           "head mask) or an array with the dimensions of `volume`.")
    }
    mask <- n4_default_mask(x, d)
    if (verbose) {
      message(sprintf("[N4] automatic head mask: %d of %d voxels", sum(mask), n))
    }
  } else {
    mask <- validate_n4_array(mask, d, "mask", fname)
    mask <- is.finite(mask) & mask != 0
  }
  if (!any(mask)) {
    stop(fname, ": the mask is empty.")
  }
  if (is.null(weight_mask)) {
    weight <- numeric(0)
  } else {
    weight <- validate_n4_array(weight_mask, d, "weight_mask", fname)
    if (any(!is.finite(weight)) || any(weight < 0)) {
      stop(fname, ": `weight_mask` must be finite and non-negative.")
    }
  }

  # ---- 3. N4 ----------------------------------------------------------------
  res <- n4_bias_field(
    x = x, dims = d, mask = mask, weight = weight, spacing = as.double(spacing),
    shrink = shrink_factor, iterations = iterations, tol = as.double(tolerance),
    spline_distance = spline_distance, order = spline_order, bins = histogram_bins,
    fwhm = as.double(bias_fwhm), noise = as.double(wiener_noise), verbose = verbose)

  bias <- exp(res$log_bias)

  # ---- 4. output --------------------------------------------------------------
  if (return_bias_field) {
    dim(bias) <- d
    if (!is.null(vox2ras)) attr(bias, "vox2ras") <- vox2ras
    return(bias)
  }

  corrected <- x / bias
  if (rescale_intensities) {
    inside <- which(mask)
    original <- x[inside]
    current <- corrected[inside]
    ok <- is.finite(original) & is.finite(current)
    range_original <- range(original[ok])
    range_current <- range(current[ok])
    if (diff(range_current) > 0 && all(is.finite(range_original))) {
      slope <- diff(range_original) / diff(range_current)
      corrected[inside] <- range_original[2] - slope * (range_current[2] - current)
    }
  }
  dim(corrected) <- d
  if (!is.null(vox2ras)) attr(corrected, "vox2ras") <- vox2ras
  corrected
}


# ---- internal helpers ---------------------------------------------------------

# Automatic head mask (mask = "auto"), following ANTsPy `get_mask(image)` with cleanup = 2:
# threshold at the image mean (and at most the image max), erode by 2, keep the
# largest component, dilate by 2, fill holes; retry with less clean-up when the
# result is empty or covers everything.
n4_default_mask <- function(x, dims, cleanup = 2L) {
  finite <- is.finite(x)
  if (!any(finite)) {
    stop("`bias_correction_n4`: `volume` has no finite voxel.")
  }
  threshold <- mean(x[finite])
  upper <- max(x[finite])
  base <- finite & x >= threshold & x <= upper
  mask <- base
  cleanup <- as.integer(cleanup)
  if (cleanup > 0L) {
    mask <- n4_mask_cleanup(base, dims, cleanup, TRUE, 50L)
    while (cleanup > 0L && (all(mask) || !any(mask))) {
      cleanup <- cleanup - 1L
      mask <- base
      if (cleanup > 0L) {
        mask <- n4_mask_cleanup(base, dims, cleanup, FALSE, 50L)
      }
    }
  }
  mask
}

validate_n4_integer <- function(v, name, minimum, fname) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < minimum || v != round(v)) {
    stop(fname, ": `", name, "` must be a single integer >= ", minimum, ".")
  }
  as.integer(v)
}

# An optional companion array (mask / weight) as a double vector matching the
# volume dimensions.
validate_n4_array <- function(a, dims, name, fname) {
  da <- dim(a)
  if (is.null(da) || length(da) < 3 || !all(da[1:3] == dims) ||
      (length(da) > 3 && !all(da[-(1:3)] == 1))) {
    stop(fname, ": `", name, "` must be an array with the same dimensions as `volume` (",
         paste(dims, collapse = " x "), ").")
  }
  as.double(a)
}
