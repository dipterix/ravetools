#' @title Prior-guided tissue segmentation with a \verb{Gaussian} mixture and an \verb{MRF}
#' @description
#' Segments a 3D intensity volume into \code{K} tissue classes with a finite
#' \verb{Gaussian} mixture model, optional spatial prior probability maps, and
#' a mean-field \verb{Markov} random field (\verb{MRF}) that favors
#' spatially coherent labels. This is a self-contained re-implementation of
#' the \verb{Atropos} algorithm of \code{'ANTs'} (\verb{Avants} and
#' colleagues, 2011) in \pkg{Rcpp}; the defaults reproduce the
#' \code{'ANTsPy'} \code{atropos}
#' call \code{i = 'Kmeans[3]'}, \code{m = '[0.2,1x1x1]'}, \code{c = '[5,0]'},
#' \code{priorweight = 0.25}.
#'
#' @param volume a 3D numeric (or integer, logical) array of intensities;
#' a \code{"vox2ras"} attribute, if present, is used when \code{vox2ras} is
#' \code{NULL}. Trailing singleton dimensions (for example a
#' \verb{NIfTI} volume read as \code{nx x ny x nz x 1}) are dropped from
#' \code{volume}, \code{mask} and \code{priors}, and the outputs are always
#' 3D. Extreme outlier intensities inside the mask (hot voxels, raw
#' \verb{CT} values) should be truncated or clipped beforehand, since the
#' \verb{k-means} seeds are spread over the intensity range (see
#' 'Initialization'). The result also depends on the absolute intensity
#' scale, because the \verb{Gaussian} density is floored at \code{1e-10} as
#' in \code{'ANTs'} (see 'Details'): intensities in very large units (class
#' standard deviations beyond roughly \code{1e7}) should be rescaled first
#' @param mask optional 3D array of the same dimensions as \code{volume};
#' non-zero (and non-\code{NA}) voxels are segmented, everything else
#' receives label \code{0} and zero posteriors. The default \code{NULL} uses
#' every voxel with a finite intensity. Voxels with non-finite intensities
#' are always excluded from the mask
#' @param priors optional list of \code{K} non-negative 3D arrays on the
#' grid of \code{volume}, the spatial prior probability of each class (for
#' example tissue probability maps warped from a template with
#' \code{apply_transform3d_volume}). They are normalized per voxel; at an
#' isolated voxel where all priors are zero a uniform prior is used, but
#' priors that are zero at every voxel inside the mask (for example
#' integer-truncated probability maps, or maps warped onto the wrong grid)
#' are an error, as is a prior that is zero everywhere inside the mask. The
#' number of classes is the number of priors, and the output classes follow
#' the order of the list. Values outside the mask are ignored
#' @param n_classes number of classes for the \verb{k-means} initialization
#' (default \code{3}, the classic \verb{CSF}, gray matter, white matter);
#' ignored when \code{priors} are given, unless it is passed explicitly and
#' conflicts with the number of priors, which is an error
#' @param prior_weight weight \eqn{w \in [0, 1]} of the spatial priors in the
#' posterior (default \code{0.25}), with the meaning of the
#' \code{priorweight} of \code{'ANTs'}; \code{0} uses the priors only to
#' initialize the model, \code{1} ignores the image and segments by the
#' priors alone (see 'Details'). Ignored without \code{priors}
#' @param mrf_beta \verb{MRF} smoothing factor \eqn{\beta \ge 0} (default
#' \code{0.2}); larger values produce smoother label maps, \code{0} disables
#' the \verb{MRF}
#' @param mrf_radius \verb{MRF} neighborhood radius in voxels, either one
#' non-negative integer for all axes or one per axis (default \code{1}, the
#' 26-connected neighborhood); non-integer values are an error and \code{0}
#' disables the \verb{MRF}. A radius beyond the grid extent is equivalent to
#' the extent along that axis and is clamped to it
#' @param iterations maximum number of \verb{EM} iterations (default
#' \code{5}); the iterations usually stop earlier, see \code{tolerance}
#' @param tolerance non-negative convergence threshold on the mean maximum
#' posterior probability over the mask (default \code{0}), with the meaning
#' of the convergence threshold of \code{'ANTs'} \verb{Atropos}: after each
#' iteration (except the first) the change of this value from the previous
#' iteration is compared with \code{tolerance}, and the iterations stop as
#' soon as the change is below it. With the default \code{0} the iterations
#' stop as soon as the value decreases (as \verb{Atropos} does with
#' \code{c = [n,0]}), otherwise they run up to \code{iterations}; the
#' posteriors and labels of the last iteration performed are returned (see
#' 'Details'). Negative values are an error
#' @param vox2ras optional 4x4 (or 3x4) matrix mapping 0-indexed voxel
#' coordinates to \verb{RAS}; only its \eqn{3\times 3} part matters, to
#' express the \verb{MRF} neighbor weights in physical distance. Default
#' \code{NULL} looks for a \code{"vox2ras"} attribute on \code{volume} and
#' otherwise uses voxel units (unit spacing)
#' @param verbose logical; print per-iteration progress (default
#' \code{FALSE})
#'
#' @details
#' \strong{Model.} For a voxel \eqn{i} inside the mask with intensity
#' \eqn{y_i}, the posterior of class \eqn{k} follows the default
#' (\code{Socrates}) posterior formulation of \verb{Atropos},
#' \deqn{P_i(k) \propto s_k(i)^{w}\,\left(N(y_i \mid \mu_k, \sigma_k^2)\,
#' M_k(i)\right)^{1-w},}{%
#' P_i(k) proportional to s_k(i)^w * (N(y_i | mu_k, sigma_k^2) * M_k(i))^(1 - w),}
#' where \eqn{w} is \code{prior_weight},
#' \deqn{s_k(i) = \frac{\pi_k\, p_k(i)}{\sum_c \pi_c\, p_c(i)}}{%
#' s_k(i) = pi_k p_k(i) / sum_c pi_c p_c(i)}
#' is the spatial prior (the per-voxel normalized prior \eqn{p_k(i)}
#' re-weighted by the mixture proportions \eqn{\pi_k}{pi_k}) and
#' \deqn{M_k(i) = \frac{\exp\left(\beta \sum_j \omega_{ij} P_j(k)\right)}{%
#' \sum_c \exp\left(\beta \sum_j \omega_{ij} P_j(c)\right)}}{%
#' M_k(i) = exp(beta * sum_j omega_ij P_j(k)) / sum_c exp(beta * sum_j omega_ij P_j(c))}
#' is the \verb{MRF} term, with \eqn{\beta}{beta} equal to \code{mrf_beta},
#' the sums over the \eqn{(2r+1)^3 - 1} neighbors \eqn{j} inside the mask, their
#' posteriors \eqn{P_j(k)} and the inverse physical distance
#' \eqn{\omega_{ij} = 1 / d_{ij}}{omega_ij = 1 / d_ij} between the two voxels.
#' As in \code{'ANTs'}, \eqn{s_k}, \eqn{M_k} and the \verb{Gaussian} density
#' are each floored at \code{1e-10} before the powers are taken, so a zero
#' prior is a very strong but not absolute constraint. Consequently, with
#' \code{prior_weight = 0} the priors only initialize the model; with
#' \code{prior_weight = 1} the image and the \verb{MRF} are ignored and the
#' segmentation is \eqn{\arg\max_k \pi_k p_k(i)}{argmax_k pi_k p_k(i)}; and
#' without \code{priors} the spatial prior is the same constant for every
#' class (\verb{Atropos} applies no prior weight then), so the posterior is
#' proportional to \eqn{N(y_i \mid \mu_k, \sigma_k^2)\, M_k(i)}{N(y_i | mu_k,
#' sigma_k^2) * M_k(i)} and the mixture proportions do not enter it. Each
#' iteration performs one synchronous E-step (every voxel reads the
#' posteriors of the previous iteration, initially the one-hot initial
#' labels) followed by an M-step that re-estimates \eqn{\mu_k}{mu_k} and
#' \eqn{\sigma_k}{sigma_k} (posterior-weighted mean and unbiased weighted
#' variance, the \verb{ITK} estimators used by \verb{Atropos}, with a small
#' variance floor relative to the intensity variance inside the mask) and
#' \eqn{\pi_k}{pi_k} (the mean posterior over the mask).
#'
#' \strong{Intensity scale.} Because the floor of \code{1e-10} is applied to
#' the (unnormalized) \verb{Gaussian} density, whose peak is
#' \eqn{1 / (\sqrt{2\pi}\,\sigma_k)}{1 / (sqrt(2 pi) sigma_k)}, the result is
#' not invariant to the intensity scale, exactly as in \verb{Atropos}: once a
#' class standard deviation exceeds roughly \code{1e7} the floor starts to
#' flatten its likelihood, and beyond roughly \code{1e9} the \verb{EM}
#' inflates the variances until every likelihood is floored and the
#' segmentation collapses to a single class. Intensities in such units
#' should be rescaled (or truncated) before segmenting; the usual \verb{MRI}
#' and \verb{CT} ranges are unaffected.
#'
#' \strong{Convergence.} After each iteration the mean maximum posterior
#' probability over the mask (the returned \code{trace}) is compared with the
#' previous iteration's value, and the iterations stop as soon as the signed
#' change is below \code{tolerance}; the first iteration never stops. This is
#' the stopping rule of \verb{Atropos}, so with the default
#' \code{tolerance = 0} the iterations stop after the first iteration in
#' which this value decreased, and the posteriors, labels and parameters
#' of that iteration are returned (there is no rollback to the previous
#' iteration). The value itself equals the posterior probability that
#' \verb{Atropos} reports only without \code{priors} or with
#' \code{prior_weight = 0}, and with \code{mrf_beta = 0}: \verb{Atropos}
#' computes its measure during its iterated-conditional-modes sweep, with the
#' neighbors' hard labels and the normalized prior \eqn{p_k(i)} itself (not
#' re-weighted by the proportions) instead of \eqn{s_k(i)}, so with weighted
#' priors or an active \verb{MRF} the two measures, and hence the iteration
#' at which the two implementations stop, can differ.
#'
#' \strong{Initialization.} Without \code{priors}, a deterministic
#' \verb{k-means} (Lloyd iterations seeded at equally spaced intensities
#' between the minimum and the maximum inside the mask, as in \code{'ANTs'})
#' clusters the masked intensities and the classes are ordered by ascending
#' mean, so with three classes label \code{1} is the darkest tissue. A
#' handful of extreme intensities stretches that range so much that a seed
#' can end up without any voxel, which is reported as an error naming the
#' intensity range: truncate or clip the intensities first (the intensity
#' truncation of \code{bias_correction_n4} does this), tighten
#' the mask, or reduce \code{n_classes}. With \code{priors}, each voxel
#' starts at the class with the largest prior (voxels where every prior is
#' zero start unlabeled), the initial class statistics are weighted by that
#' prior value and the initial proportions are the label fractions; the
#' class order is the order of the list.
#'
#' \strong{Relation to \code{'ANTs'} \verb{Atropos}.} The posterior formula
#' above, the \verb{Gaussian} likelihood, the \verb{k-means} initialization,
#' the inverse-distance neighborhood weights, the probability floor, the
#' parameter estimators and the stopping rule are those of \verb{Atropos}
#' with its defaults (the \code{Socrates} formulation with mixture
#' proportions and no annealing), so \code{prior_weight} is interchangeable
#' with its \code{priorweight}. The remaining difference is the \verb{MRF}
#' term: \verb{Atropos} plugs the neighbors' hard labels (updated by an
#' iterated-conditional-modes sweep) into \eqn{M_k(i)}, whereas here the
#' neighbors' posteriors are used (a mean-field update), so the two can
#' differ along tissue boundaries when \code{mrf_beta > 0}; with
#' \code{mrf_beta = 0} the same model is computed, although the convergence
#' measure (see above) and therefore the number of iterations can still
#' differ when \code{prior_weight > 0}. The posteriors returned here are the
#' ones the labels were taken from, whereas \verb{Atropos} re-evaluates the
#' probability images it writes out with the parameters updated in its last
#' iteration. Two safeguards also differ: a class that receives no voxel at
#' initialization falls back to prior-weighted statistics instead of a zero
#' likelihood, and priors that are zero everywhere inside the mask are an
#' error.
#'
#' \strong{Memory.} All computation is restricted to the mask; the working
#' memory is about \eqn{2 \times 8 K N_{mask}}{2 * 8 * K * N_mask} bytes for
#' the two posterior buffers (plus \eqn{8 K N_{mask}}{8 * K * N_mask} when
#' \code{prior_weight > 0}) and one integer per voxel of the full grid.
#'
#' @returns A list with
#' \describe{
#' \item{\code{segmentation}}{integer 3D array, \code{0} outside the mask
#' and \code{1..K} inside, the \verb{argmax} of the posteriors}
#' \item{\code{posteriors}}{list of \code{K} double 3D arrays with the
#' posterior probability of each class (summing to one inside the mask,
#' \code{0} outside)}
#' \item{\code{means}, \code{sds}, \code{proportions}}{the final class
#' intensity means, standard deviations and mixture proportions, estimated
#' from the returned posteriors}
#' \item{\code{trace}}{the mean maximum posterior probability over the mask
#' after each iteration performed (one entry per iteration, so its length is
#' the number of iterations run; see \code{tolerance})}
#' }
#' The \code{segmentation} and each posterior carry the \code{"vox2ras"}
#' attribute when one is known.
#'
#' @references
#' \verb{Avants} BB, \verb{Tustison} NJ, Wu J, Cook PA, Gee \verb{JC} (2011).
#' An open source multivariate framework for n-tissue segmentation with
#' evaluation on public data. \emph{\verb{Neuroinformatics}}, 9(4), 381-400.
#' \doi{10.1007/s12021-011-9109-y}
#'
#' \verb{Zhang} Y, Brady M, Smith S (2001). Segmentation of brain MR images
#' through a hidden Markov random field model and the
#' expectation-maximization algorithm. \emph{IEEE Transactions on Medical
#' Imaging}, 20(1), 45-57. \doi{10.1109/42.906424}
#'
#' \verb{Dempster} AP, Laird NM, Rubin DB (1977). Maximum likelihood from
#' incomplete data via the EM algorithm. \emph{Journal of the Royal
#' Statistical Society, Series B}, 39(1), 1-38.
#' \doi{10.1111/j.2517-6161.1977.tb01600.x}
#'
#' @seealso \code{\link{register_volume3d}} to bring template priors into the
#' subject space
#'
#' @examples
#'
#' # A toy phantom: nested spheres with three intensity classes
#' nd <- 32
#' ctr <- (nd - 1) / 2
#' idx <- arrayInd(seq_len(nd^3), rep(nd, 3)) - 1
#' r <- sqrt(rowSums((idx - ctr)^2))
#' truth <- integer(nd^3)
#' truth[r <= 14] <- 1L
#' truth[r <= 10] <- 2L
#' truth[r <= 6] <- 3L
#' dim(truth) <- rep(nd, 3)
#' mask <- truth > 0
#' set.seed(1)
#' volume <- array(0, dim(truth))
#' volume[mask] <- c(30, 80, 130)[truth[mask]] + rnorm(sum(mask), sd = 12)
#'
#' # k-means initialization, classes ordered by intensity
#' res <- segment_volume_tissue_gmm(volume, mask = mask)
#' res$means
#' table(truth = truth[mask], segmented = res$segmentation[mask])
#'
#' # the same with spatial priors (here a softened version of the truth, in a
#' # shuffled order) that fix the class order
#' priors <- lapply(c(3, 1, 2), function(k) {
#'   p <- (truth == k) * 1
#'   p <- (p + 0.1) / 1.3   # a soft, not quite informative prior
#'   p
#' })
#' res2 <- segment_volume_tissue_gmm(volume, mask = mask, priors = priors,
#'                                   prior_weight = 0.25)
#' res2$means       # follows the prior order: bright, dark, medium
#'
#' @export
segment_volume_tissue_gmm <- function(
    volume, mask = NULL, priors = NULL, n_classes = 3L,
    prior_weight = 0.25, mrf_beta = 0.2, mrf_radius = 1L, iterations = 5L,
    tolerance = 0, vox2ras = NULL, verbose = FALSE) {

  # `n_classes = NULL` counts as not given (the priors then define K)
  n_classes_given <- !missing(n_classes) && !is.null(n_classes)
  fname <- "`segment_volume_tissue_gmm`"

  if (is.null(vox2ras)) vox2ras <- attr(volume, "vox2ras")
  volume <- gmm_as_volume(volume, "volume", fname)
  d <- dim(volume)

  # Mask: non-zero and non-NA voxels, always restricted to finite intensities
  if (is.null(mask)) {
    mask <- is.finite(volume)
  } else {
    mask <- gmm_as_volume(mask, "mask", fname)
    if (!identical(dim(mask), d)) {
      stop(sprintf("%s: `mask` must have the same dimensions as `volume` (%s).",
                   fname, paste(d, collapse = " x ")))
    }
    mask <- (mask != 0) & !is.na(mask) & is.finite(volume)
  }
  dim(mask) <- NULL
  n_mask <- sum(mask)
  if (n_mask == 0L) {
    stop(sprintf("%s: the mask contains no voxel with a finite intensity.", fname))
  }

  # `n_classes` must be a single integer in [2, .Machine$integer.max] whenever
  # it is used: always without priors, and with priors only when given
  if (n_classes_given || is.null(priors)) {
    if (!is.numeric(n_classes) || length(n_classes) != 1L || !is.finite(n_classes) ||
        n_classes != round(n_classes) || n_classes < 2 ||
        n_classes > .Machine$integer.max) {
      stop(sprintf("%s: `n_classes` must be a single integer >= 2 (at most %d).",
                   fname, .Machine$integer.max))
    }
  }

  # Priors define the number and order of the classes
  if (!is.null(priors)) {
    if (!is.list(priors) || length(priors) < 2L) {
      stop(sprintf("%s: `priors` must be a list of at least two 3D arrays.", fname))
    }
    k_classes <- length(priors)
    if (n_classes_given && n_classes != k_classes) {
      stop(sprintf(
        "%s: `n_classes` (%d) conflicts with the number of priors (%d); omit `n_classes` when `priors` are given.",
        fname, as.integer(n_classes), k_classes))
    }
    priors <- lapply(seq_len(k_classes), function(k) {
      p <- gmm_as_volume(priors[[k]], sprintf("priors[[%d]]", k), fname)
      if (!identical(dim(p), d)) {
        stop(sprintf("%s: `priors[[%d]]` must have the same dimensions as `volume` (%s).",
                     fname, k, paste(d, collapse = " x ")))
      }
      dim(p) <- NULL
      p
    })
  } else {
    k_classes <- as.integer(n_classes)
    priors <- list()
  }
  if (n_mask < k_classes) {
    stop(sprintf("%s: the mask has fewer voxels (%d) than classes (%d).",
                 fname, n_mask, k_classes))
  }

  # Scalar parameters: numeric (not character or logical), finite, in range
  if (!is.numeric(prior_weight) || length(prior_weight) != 1L ||
      !is.finite(prior_weight) || prior_weight < 0 || prior_weight > 1) {
    stop(sprintf("%s: `prior_weight` must be a single number in [0, 1].", fname))
  }
  prior_weight <- as.double(prior_weight)
  if (!is.numeric(mrf_beta) || length(mrf_beta) != 1L || !is.finite(mrf_beta) ||
      mrf_beta < 0) {
    stop(sprintf("%s: `mrf_beta` must be a single non-negative number.", fname))
  }
  mrf_beta <- as.double(mrf_beta)
  if (!is.numeric(mrf_radius) || !(length(mrf_radius) %in% c(1L, 3L)) ||
      any(!is.finite(mrf_radius)) || any(mrf_radius < 0) ||
      any(mrf_radius != round(mrf_radius))) {
    stop(sprintf("%s: `mrf_radius` must be one or three non-negative integers.", fname))
  }
  if (length(mrf_radius) == 1L) mrf_radius <- rep(mrf_radius, 3L)
  # radii beyond the grid extent are equivalent to the extent (clamped in C++)
  mrf_radius <- as.integer(pmin(mrf_radius, .Machine$integer.max))
  if (!is.numeric(iterations) || length(iterations) != 1L || !is.finite(iterations) ||
      iterations != round(iterations) || iterations < 1 ||
      iterations > .Machine$integer.max) {
    stop(sprintf("%s: `iterations` must be a single integer >= 1.", fname))
  }
  iterations <- as.integer(iterations)
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) ||
      tolerance < 0) {
    stop(sprintf("%s: `tolerance` must be a single non-negative number.", fname))
  }
  tolerance <- as.double(tolerance)

  # Physical neighbor distances come from the 3x3 part of vox2ras
  if (!is.null(vox2ras)) {
    vox2ras <- as.matrix(vox2ras)
    if (nrow(vox2ras) == 3L && ncol(vox2ras) == 4L) vox2ras <- rbind(vox2ras, c(0, 0, 0, 1))
    if (!is.numeric(vox2ras) || !all(dim(vox2ras) == c(4L, 4L)) || anyNA(vox2ras) ||
        !all(is.finite(vox2ras))) {
      stop(sprintf("%s: `vox2ras` must be a finite 4x4 (or 3x4) numeric matrix.", fname))
    }
    storage.mode(vox2ras) <- "double"
    direction <- vox2ras[1:3, 1:3]
    if (abs(det(direction)) <= 0) {
      stop(sprintf("%s: `vox2ras` is singular.", fname))
    }
  } else {
    direction <- diag(3)
  }

  # Errors raised by the native kernel (data problems found while it runs, such
  # as constant intensities or an empty class at initialization) are re-raised
  # with the same attribution as the checks above
  user_call <- sys.call()
  res <- tryCatch(
    segment_gmm_mrf_cpp(
      volume = as.double(volume), dims = d, mask = mask, priors = priors,
      n_classes = as.integer(k_classes), prior_weight = prior_weight,
      mrf_beta = mrf_beta, mrf_radius = mrf_radius, direction = direction,
      iterations = iterations, tolerance = tolerance, verbose = isTRUE(verbose)),
    error = function(e) {
      stop(simpleError(sprintf("%s: %s", fname, conditionMessage(e)), call = user_call))
    }
  )

  segmentation <- res$segmentation
  dim(segmentation) <- d
  posteriors <- lapply(res$posteriors, function(p) {
    dim(p) <- d
    if (!is.null(vox2ras)) attr(p, "vox2ras") <- vox2ras
    p
  })
  if (!is.null(vox2ras)) attr(segmentation, "vox2ras") <- vox2ras

  list(
    segmentation = segmentation,
    posteriors = posteriors,
    means = res$means,
    sds = res$sds,
    proportions = res$proportions,
    trace = res$trace
  )
}

# Coerce a numeric/integer/logical array to a plain double 3D array
gmm_as_volume <- function(x, name, fname) {
  d <- dim(x)
  if (is.null(d) || length(d) < 3L || !(is.numeric(x) || is.logical(x))) {
    stop(sprintf("%s: `%s` must be a numeric 3D array.", fname, name))
  }
  if (length(d) > 3L) {
    if (!all(d[-(1:3)] == 1L)) {
      stop(sprintf("%s: `%s` must be a 3D array (extra dimensions must be 1).", fname, name))
    }
    d <- d[1:3]
  }
  if (any(d < 1L)) {
    stop(sprintf("%s: `%s` has an empty dimension.", fname, name))
  }
  x <- as.double(x)
  dim(x) <- d
  x
}
