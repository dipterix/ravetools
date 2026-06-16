#' @title Native 3D volume registration (\code{'rigid'}, \code{'affine'}, or \code{'SyN'})
#' @description
#' Self-contained image registration for 3D volumes, implemented purely in
#' \pkg{RcppEigen} (no other external registration library). It mirrors the
#' core behavior of \pkg{'ANTs'} \code{antsRegistration}: a multi-resolution,
#' physical-shift scaled gradient-descent optimizer driving a similarity metric,
#' working entirely in the anatomical \verb{RAS} (right-anterior-superior)
#' space. Each volume carries its own \code{vox2ras} (0-indexed voxel index to
#' anatomical \verb{RAS}) \eqn{4\times 4} transform, so volumes with
#' different sampling, orientation, or field of view are aligned correctly.
#'
#' @param source the moving volume to be aligned, a 3D array (for example a
#' \code{'CT'}); the result transform maps \code{target} into this image's space.
#' May also be a \code{list} of already co-registered 3D arrays (multiple
#' modalities such as \code{'T1'}, \code{'T2'}, an atlas constraint); all source
#' channels must share the same dimensions and \code{vox2ras}
#' @param target the fixed/reference volume to align to, a 3D array (for example
#' a \code{'MRI'}), or a \code{list} of co-registered arrays matching
#' \code{source} channel-for-channel
#' @param weights optional numeric weights, one per source/target pair,
#' controlling each channel's contribution to the deformable cost; default is
#' equal weighting. Weights are normalized internally to sum to 1. Only the
#' deformable (\code{"syn"}/\code{"syn_only"}) stage is multivariate; the linear
#' stage always uses the first (primary) pair
#' @param source_vox2ras,target_vox2ras 4x4 (or 3x4) matrices mapping the
#' 0-indexed voxel coordinate (column-row-slice, \code{'C'}-style starting from
#' 0, complying with \code{'NIfTI'}) to the \verb{RAS} coordinate system; if
#' \code{NULL}, the function looks for a \code{"vox2ras"} attribute on the array
#' @param source_mask,target_mask optional 3D mask arrays restricting where the
#' metric is evaluated; default \code{NULL} (no mask, evaluate everywhere). A
#' \code{target_mask} (on the target/fixed grid, same dimensions as
#' \code{target}) limits which target voxels drive the registration; a
#' \code{source_mask} (on the source/moving grid) drops samples that map outside
#' it. Non-zero voxels are included. Besides focusing the alignment on a region
#' of interest, masks speed things up by skipping background: the linear stage
#' samples only inside the mask, and the deformable stage skips warping voxels
#' outside the (dilated) mask. One mask per grid, shared across channels
#' @param source_points,target_points optional \code{N x 3} matrices of
#' corresponding landmark coordinates that add a surface/landmark term to the
#' \strong{deformable} (\code{"syn"}/\code{"syn_only"}) stage; default
#' \code{NULL} (no term). Row \code{i} of \code{target_points} (in the
#' target/fixed \verb{RAS}) and row \code{i} of \code{source_points} (in the
#' source/moving \verb{RAS}) must be the \emph{same} anatomical location, for
#' example corresponding cortical-surface vertices from a \code{'FreeSurfer'}
#' spherical registration. The term pulls the warp of each target point onto its
#' source correspondent, recovering cortical folding (\verb{gyrification})
#' that the intensity metric blurs over, while the image metric still drives
#' deep-brain \verb{subcortical} structures. Must be supplied together with
#' equal row counts. \strong{Points must be in the same \verb{RAS} frame as
#' \code{source_vox2ras}/\code{target_vox2ras}} (note \code{'FreeSurfer'}
#' surfaces use \verb{surface/tkr} \verb{RAS}, which differs from scanner
#' \verb{RAS} by \code{c_ras})
#' @param points_weight relative weight of the landmark term against the image
#' metric in the deformable stage; default \code{0.5}. Larger values follow the
#' landmarks more closely. The sparse landmark force is attenuated by
#' \code{syn_sigma} smoothing, so this typically needs tuning for a given point
#' count and spacing
#' @param type type of transform to estimate; one of \code{'rigid'} (6 degrees
#' of freedom), \code{'affine'} (12), \code{"syn"} (\verb{affine} followed
#' by a \verb{SDR} - symmetric \verb{diffeomorphic} deformation), or
#' \code{'syn_only'} (\verb{deformable} stage only — no \verb{affine} is
#' estimated; \code{init_transform} is used directly as the starting
#' \verb{affine}, useful when you already have a good linear alignment)
#' @param metric similarity metric: \code{"mattes"} (Mattes mutual information,
#' the default, best for cross-modal such as \code{'CT'}-\code{'MRI'}),
#' \code{"cc"} (normalized cross-correlation, for same-modality), or
#' \code{"meansquares"} (mean squared intensity difference). With multiple
#' channels, supply either a single metric (used for every pair) or a vector
#' with one metric per source/target pair
#' @param shrink_factors integer down-sampling factors, one per resolution level
#' (coarsest first); default \code{c(4, 2, 1)}
#' @param smoothing_sigmas Gaussian smoothing applied at each level, in voxels,
#' same length as \code{shrink_factors}; default \code{c(2, 1, 0)}
#' @param iterations maximum optimizer iterations per level; default
#' \code{c(1000, 500, 250)}
#' @param sampling_rate fraction of fixed voxels sampled to evaluate the metric
#' (speeds up large volumes); default \code{0.2}
#' @param interpolation output interpolation used when warping each modality
#' onto the target grid: \code{'trilinear'} (default), \code{'nearest'} (keeps
#' label or segmentation values intact), or \code{"bspline"} (cubic
#' \verb{Catmull-Rom}, higher quality, slower). Like \code{metric}, supply a
#' single value (applied to every channel) or one per source/target pair. This
#' affects only the returned warped images, never the optimization's internal
#' sampling (always \code{'trilinear'}, so every channel still produces
#' smooth gradients)
#' @param number_of_bins number of histogram bins for the \code{"mattes"}
#' metric; default \code{32}
#' @param seed random seed for the voxel sampler, for reproducibility
#' @param init_transform optional 4x4 initial \verb{RAS}-to-\verb{RAS}
#' transform (fixed to moving) to start from
#' @param syn_iterations,syn_sigma deformable stage controls (only used when
#' \code{type = "syn"} or \code{"syn_only"}): per-level iteration counts and the
#' Gaussian regularization sigma (in voxels) applied to the update field. Both
#' are recycled to the number of levels, so \code{syn_sigma} may be a vector to
#' vary the regularization per stage (e.g. \code{c(3, 3, 1)} to relax it at the
#' finest level for sharper detail)
#' @param verbose logical; if \code{TRUE} (default) print per-level and
#' per-iteration progress to the console, including the current stage, shrink
#' factor, smoothing sigma, cost metric, and step size (linear) or maximum
#' displacement (deformable); useful to monitor convergence on large volumes
#' @returns A list with:
#' \describe{
#' \item{\code{transform}}{the estimated 4x4 \verb{RAS}-to-\verb{RAS} linear
#' transform mapping \code{target} (fixed) coordinates to \code{source} (moving)
#' coordinates}
#' \item{\code{image}}{the (primary) \code{source} resampled onto the
#' \code{target} grid}
#' \item{\code{images}}{a list with every source channel resampled onto the
#' \code{target} grid using its own \code{interpolation}; \code{image} is the
#' first element. For single-image input this is a length-1 list}
#' \item{\code{forward_field},\code{inverse_field}}{(only for \code{"syn"}) the
#' deformation fields}
#' \item{\code{metric_trace}}{the metric value across optimizer iterations}
#' \item{\code{type},\code{metric}}{echoes of the inputs}
#' }
#' @seealso \code{\link{apply_transform3d}}, \code{\link{resample_3d_volume}}
#' @examples
#' \donttest{
#'
#' # --- synthetic same-modality example -------------------------------------
#' nd <- c(50, 50, 50)
#' vox2ras <- diag(4); vox2ras[1:3, 4] <- -25
#' blob <- function(cx, cy, cz, s = 6) {
#'   g <- expand.grid(x = 0:(nd[1]-1), y = 0:(nd[2]-1), z = 0:(nd[3]-1))
#'   array(exp(-((g$x-cx)^2 + (g$y-cy)^2 + (g$z-cz)^2) / (2*s^2)), nd)
#' }
#' target <- blob(25, 25, 25)
#' source <- blob(28, 23, 26)            # shifted by a known (3, -2, 1) mm
#'
#' res <- register_volume3d(
#'   source, target,
#'   source_vox2ras = vox2ras, target_vox2ras = vox2ras,
#'   type = "rigid", metric = "cc"
#' )
#' res$transform[1:3, 4]                 # ~ c(3, -2, 1)
#'
#'
#' # --- multimodal registration (several co-registered channels) ------------
#' # Two aligned modalities (e.g. a T1 and a T2) jointly drive the deformable
#' # stage. Each pair may use its own metric, and `weights` set their relative
#' # influence (normalized internally to sum to 1). All source channels must
#' # share a grid; likewise all target channels.
#' t1_target <- blob(25, 25, 25, s = 6)
#' t2_target <- blob(25, 25, 25, s = 9)  # same anatomy, different contrast
#' t1_source <- blob(27, 24, 26, s = 6)  # both channels share the same warp
#' t2_source <- blob(27, 24, 26, s = 9)
#'
#' res_mm <- register_volume3d(
#'   source = list(t1_source, t2_source),
#'   target = list(t1_target, t2_target),
#'   source_vox2ras = vox2ras, target_vox2ras = vox2ras,
#'   weights = c(2, 1),                  # T1 counts twice as much as T2
#'   metric = c("mattes", "cc"),         # one metric per channel
#'   type = "syn"
#' )
#' res_mm$transform[1:3, 4]
#'
#' }
#' @export
register_volume3d <- function(
    source, target,
    source_vox2ras = NULL, target_vox2ras = NULL,
    source_mask = NULL, target_mask = NULL,
    source_points = NULL, target_points = NULL, points_weight = 0.5,
    weights = NULL,
    type = c("rigid", "affine", "syn", "syn_only"),
    metric = "mattes",
    shrink_factors = c(4, 2, 1),
    smoothing_sigmas = c(2, 1, 0),
    iterations = c(1000, 500, 250),
    sampling_rate = 0.2,
    interpolation = "trilinear",
    number_of_bins = 32L,
    seed = 1L,
    init_transform = NULL,
    syn_iterations = c(40, 20, 0),
    syn_sigma = 3,
    verbose = TRUE) {

  type <- match.arg(type)

  # Accept either a single array or a list of co-registered modalities.
  if (!is.list(source)) source <- list(source)
  if (!is.list(target)) target <- list(target)
  n_pairs <- length(source)
  if (length(target) != n_pairs) {
    stop("`register_volume3d`: `source` and `target` must have the same number of channels.")
  }

  # Per-channel metric. A single metric (the default "mattes") is recycled to
  # every channel; otherwise supply one metric per pair. Each entry is validated
  # against the available options, so three channels left at the default all
  # use "mattes" rather than being silently read as c("mattes","cc","meansquares").
  metric_choices <- c("mattes", "cc", "meansquares")
  metric <- as.character(metric)
  if (length(metric) == 1L) metric <- rep(metric, n_pairs)
  if (length(metric) != n_pairs) {
    stop("`register_volume3d`: `metric` must have length 1 or length(source).")
  }
  if (!all(metric %in% metric_choices)) {
    stop("`register_volume3d`: `metric` must be one of ",
         paste(metric_choices, collapse = ", "), ".")
  }

  # Per-channel output interpolation, validated like `metric`. Controls how each
  # warped modality is resampled onto the target grid, e.g. "nearest" to keep a
  # segmentation's label values intact. The optimizer's internal metric sampling
  # is always trilinear, so this does not change what drives the registration.
  interp_choices <- c("trilinear", "nearest", "bspline")
  interpolation <- as.character(interpolation)
  if (length(interpolation) == 1L) interpolation <- rep(interpolation, n_pairs)
  if (length(interpolation) != n_pairs) {
    stop("`register_volume3d`: `interpolation` must have length 1 or length(source).")
  }
  if (!all(interpolation %in% interp_choices)) {
    stop("`register_volume3d`: `interpolation` must be one of ",
         paste(interp_choices, collapse = ", "), ".")
  }
  interp_codes <- vapply(interpolation,
    function(x) switch(x, nearest = 0L, bspline = 2L, 1L), integer(1))

  # Weights: default equal, normalized to sum to 1.
  if (is.null(weights)) weights <- rep(1.0, n_pairs)
  weights <- as.double(weights)
  if (length(weights) != n_pairs || anyNA(weights) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("`register_volume3d`: `weights` must be non-negative, length(source), and not all zero.")
  }
  weights <- weights / sum(weights)

  # A single vox2ras applies to every channel in each group (co-registered).
  source_vox2ras <- resolve_vox2ras(source[[1]], source_vox2ras, "source")
  target_vox2ras <- resolve_vox2ras(target[[1]], target_vox2ras, "target")

  # Convert each channel; all channels in a group must share the same grid.
  src_list <- lapply(source, as_double_volume, "source")
  tgt_list <- lapply(target, as_double_volume, "target")
  src_dim <- dim(src_list[[1]])
  tgt_dim <- dim(tgt_list[[1]])
  for (k in seq_len(n_pairs)) {
    if (!identical(dim(src_list[[k]]), src_dim)) {
      stop("`register_volume3d`: all `source` channels must share the same dimensions.")
    }
    if (!identical(dim(tgt_list[[k]]), tgt_dim)) {
      stop("`register_volume3d`: all `target` channels must share the same dimensions.")
    }
  }

  # Optional masks (NULL = none, unchanged behavior). A mask restricts where the
  # metric is evaluated and skips background voxels: `target_mask` (fixed/target
  # grid) limits which target voxels drive the registration; `source_mask`
  # (moving/source grid) drops samples that map outside it. Nonzero = include.
  # One mask per grid (shared across channels, which are co-registered).
  target_mask <- validate_mask(target_mask, tgt_dim, "target_mask")
  source_mask <- validate_mask(source_mask, src_dim, "source_mask")

  # Optional cortical landmark correspondences (only used by the deformable
  # stage). `target_points` (fixed/target RAS) and `source_points` (moving/source
  # RAS) are N x 3 matrices whose row i is the SAME cortical vertex in each space
  # (e.g. corresponding vertices from FreeSurfer spherical registration). They
  # add a weighted term that pulls the warp of each target point onto its source
  # correspondent, recovering gyrification the intensity metric blurs over.
  # IMPORTANT: points must be in the same RAS frame as the vox2ras matrices.
  source_points <- validate_points(source_points, "source_points")
  target_points <- validate_points(target_points, "target_points")
  if (is.null(source_points) != is.null(target_points)) {
    stop("`register_volume3d`: `source_points` and `target_points` must be supplied together.")
  }
  if (!is.null(source_points) && nrow(source_points) != nrow(target_points)) {
    stop("`register_volume3d`: `source_points` and `target_points` must have the same number of rows (corresponding vertices).")
  }
  points_weight <- as.double(points_weight)[[1L]]

  # The linear (rigid/affine) stage is driven by the primary pair only; extra
  # channels contribute only to the deformable stage (ANTs convention).
  if (n_pairs > 1 && type %in% c("rigid", "affine")) {
    warning("`register_volume3d`: multiple channels only affect the deformable ",
            "(syn / syn_only) stage; using the first pair for type = '", type, "'.")
  }
  src <- src_list[[1]]
  tgt <- tgt_list[[1]]
  metric1 <- metric[[1]]

  shrink_factors <- as.integer(shrink_factors)
  nlev <- length(shrink_factors)
  smoothing_sigmas <- recycle(as.double(smoothing_sigmas), nlev)
  iterations <- recycle(as.integer(iterations), nlev)

  # Default initialization: align the centers of mass of the primary pair in
  # RAS. Without this, volumes whose scanner origins differ (e.g. CT vs MRI from
  # different sessions) start with no overlap and the optimizer cannot move.
  if (is.null(init_transform)) {
    com_fixed <- center_of_mass(tgt, target_vox2ras)
    com_moving <- center_of_mass(src, source_vox2ras)
    init_transform <- diag(4)
    init_transform[1:3, 4] <- com_moving - com_fixed
  } else {
    init_transform <- as.matrix(init_transform)
  }

  run_linear <- function(linear_type, init) {
    register_linear_cpp(
      fixed = as.double(tgt), fixedDim = tgt_dim, fixedVox2Ras = target_vox2ras,
      moving = as.double(src), movingDim = src_dim, movingVox2Ras = source_vox2ras,
      type = linear_type, metric = metric1,
      shrinkFactors = shrink_factors, smoothingSigmas = smoothing_sigmas,
      iterations = iterations, samplingRate = as.double(sampling_rate),
      learningRate = 0, numberOfBins = as.integer(number_of_bins),
      seed = as.integer(seed), initTransform = init,
      fixedMask = target_mask, movingMask = source_mask,
      verbose = isTRUE(verbose))
  }

  syn_iters <- recycle(as.integer(syn_iterations), nlev)
  syn_sigma <- recycle(as.double(syn_sigma), nlev)   # per-level flow regularization

  # syn_only: skip all linear stages, use init_transform as the affine init
  if (type == "syn_only") {
    transform <- init_transform
    syn <- register_syn(
      tgt_list, tgt_dim, target_vox2ras, src_list, src_dim, source_vox2ras,
      metrics = metric, weights = weights, init_transform = transform,
      shrink_factors = shrink_factors, smoothing_sigmas = smoothing_sigmas,
      iterations = syn_iters, interp_codes = interp_codes,
      field_sigma = syn_sigma, seed = as.integer(seed),
      fixed_mask = target_mask, moving_mask = source_mask,
      fixed_points = target_points, moving_points = source_points,
      points_weight = points_weight,
      verbose = isTRUE(verbose))
    return(finalize_registration(list(
      transform = transform,
      forward_field = syn$forward_field,
      inverse_field = syn$inverse_field,
      image = syn$image,
      images = syn$images,
      metric_trace = syn$metric_trace,
      type = type,
      metric = metric
    ), source_vox2ras, target_vox2ras, tgt_dim))
  }

  # Staged initialization (ANTs-style): rigid -> affine -> deformable
  init <- init_transform
  trace <- numeric(0)

  if (type %in% c("affine", "syn")) {
    rigid <- run_linear("rigid", init)
    init <- rigid$transform
    trace <- c(trace, rigid$metric_trace)
    linear <- run_linear("affine", init)
  } else {
    linear <- run_linear("rigid", init)
  }
  trace <- c(trace, linear$metric_trace)
  transform <- linear$transform

  result <- list(
    transform = transform,
    type = type,
    metric = metric,
    metric_trace = trace
  )

  if (type == "syn") {
    syn <- register_syn(
      tgt_list, tgt_dim, target_vox2ras, src_list, src_dim, source_vox2ras,
      metrics = metric, weights = weights, init_transform = transform,
      shrink_factors = shrink_factors, smoothing_sigmas = smoothing_sigmas,
      iterations = syn_iters, interp_codes = interp_codes,
      field_sigma = syn_sigma, seed = as.integer(seed),
      fixed_mask = target_mask, moving_mask = source_mask,
      fixed_points = target_points, moving_points = source_points,
      points_weight = points_weight,
      verbose = isTRUE(verbose))
    result$forward_field <- syn$forward_field
    result$inverse_field <- syn$inverse_field
    result$metric_trace <- c(trace, syn$metric_trace)
    result$image <- syn$image
    result$images <- syn$images
    return(finalize_registration(result, source_vox2ras, target_vox2ras, tgt_dim))
  }

  # resample each moving channel onto the fixed grid using the linear transform,
  # each with its own output interpolation
  images <- lapply(seq_len(n_pairs), function(k) {
    apply_transform3d(
      src_list[[k]], source_vox2ras, transform,
      reference_dim = tgt_dim, reference_vox2ras = target_vox2ras,
      interpolation = interpolation[k])
  })
  result$image <- images[[1]]
  result$images <- images

  finalize_registration(result, source_vox2ras, target_vox2ras, tgt_dim)
}

# Tag a registration result with its class, attach the fixed-grid vox2ras to the
# deformation fields (so they are self-describing for save_registration), and
# record the grid geometry needed to re-apply or export the transform.
finalize_registration <- function(result, source_vox2ras, target_vox2ras, target_dim) {
  if (!is.null(result$forward_field)) {
    attr(result$forward_field, "vox2ras") <- target_vox2ras
  }
  if (!is.null(result$inverse_field)) {
    attr(result$inverse_field, "vox2ras") <- target_vox2ras
  }
  result$geometry <- list(
    source_vox2ras = source_vox2ras,
    target_vox2ras = target_vox2ras,
    target_dim = target_dim
  )
  class(result) <- "ravetools_register_volume3d"
  result
}


#' @title Apply a linear \verb{RAS} transform to resample a 3D volume
#' @description
#' Warps a moving volume onto a reference grid given a 4x4
#' \verb{RAS}-to-\verb{RAS} transform (such as the \code{transform} returned
#' by \code{\link{register_volume3d}}). The transform maps reference (fixed)
#' \verb{RAS} coordinates to moving \verb{RAS} coordinates.
#' @param volume moving 3D array to resample
#' @param vox2ras the moving volume's voxel-to-\verb{RAS} 4x4 transform
#' @param transform 4x4 \verb{RAS}-to-\verb{RAS} transform (fixed to moving)
#' @param reference_dim output dimension (the fixed grid); defaults to
#' \code{dim(volume)}
#' @param reference_vox2ras the fixed grid's voxel-to-\verb{RAS} transform;
#' defaults to \code{vox2ras}
#' @param interpolation \code{'trilinear'} (default), \code{'nearest'}, or
#' \code{'bspline'} (cubic \verb{Catmull-Rom})
#' @param na_fill value for out-of-bounds voxels; default \code{0}
#' @returns The resampled volume on the reference grid, with a \code{'vox2ras'}
#' attribute equal to \code{reference_vox2ras}.
#' @seealso \code{\link{register_volume3d}}
#' @export
apply_transform3d <- function(
    volume, vox2ras, transform,
    reference_dim = dim(volume), reference_vox2ras = vox2ras,
    interpolation = c("trilinear", "nearest", "bspline"), na_fill = 0) {

  interpolation <- match.arg(interpolation)
  interp_code <- switch(interpolation, trilinear = 1L, bspline = 2L, 0L)

  volume <- as_double_volume(volume, "volume")
  vox2ras <- as.matrix(vox2ras)
  reference_vox2ras <- as.matrix(reference_vox2ras)
  transform <- as.matrix(transform)

  # for a fixed voxel f: moving voxel = solve(vox2ras) %*% transform %*%
  #   reference_vox2ras %*% f  =>  newVoxToWorld = transform %*% reference_vox2ras
  new_vox_to_world <- transform %*% reference_vox2ras

  storage.mode(na_fill) <- "double"
  re <- resample3D(
    arrayDim = as.integer(reference_dim[1:3]),
    fromArray = volume,
    newVoxToWorldTransposed = t(new_vox_to_world),
    oldVoxToWorldTransposed = t(vox2ras),
    na = na_fill,
    interpolation = interp_code)

  out <- re[[1]]
  attr(out, "vox2ras") <- reference_vox2ras
  out
}


# ---- internal helpers -------------------------------------------------------

resolve_vox2ras <- function(x, v, name) {
  if (is.null(v)) v <- attr(x, "vox2ras")
  if (is.null(v)) {
    stop(sprintf(
      "`register_volume3d`: %s_vox2ras is missing. Provide a 4x4 voxel-to-RAS matrix (or set a 'vox2ras' attribute on `%s`).",
      name, name))
  }
  v <- as.matrix(v)
  if (nrow(v) == 3 && ncol(v) == 4) v <- rbind(v, c(0, 0, 0, 1))
  if (!all(dim(v) == c(4, 4))) {
    stop(sprintf("`register_volume3d`: %s_vox2ras must be a 4x4 (or 3x4) matrix.", name))
  }
  v
}

as_double_volume <- function(x, name) {
  d <- dim(x)
  if (length(d) < 3 || any(d[1:3] <= 1)) {
    stop(sprintf("`register_volume3d`: `%s` must be a 3D volume.", name))
  }
  if (length(d) > 3 && !all(d[-(1:3)] == 1)) {
    stop(sprintf("`register_volume3d`: `%s` has more than 3 non-singleton dimensions.", name))
  }
  storage.mode(x) <- "double"
  dim(x) <- d[1:3]
  x
}

# Validate an optional mask (NULL passes through) and return it as a flat double
# vector matching `ref_dim` (the corresponding image grid). Nonzero = include.
validate_mask <- function(m, ref_dim, name) {
  if (is.null(m)) return(NULL)
  m <- as_double_volume(m, name)
  if (!identical(dim(m), ref_dim)) {
    stop(sprintf("`register_volume3d`: `%s` must match its image dimensions (%s).",
                 name, paste(ref_dim, collapse = " x ")))
  }
  as.double(m)
}

# Validate an optional landmark point set (NULL passes through) and return it as
# an N x 3 double matrix of RAS coordinates.
validate_points <- function(p, name) {
  if (is.null(p)) return(NULL)
  p <- as.matrix(p)
  storage.mode(p) <- "double"
  if (ncol(p) != 3L) {
    stop(sprintf("`register_volume3d`: `%s` must be an N x 3 matrix of RAS coordinates.", name))
  }
  p
}

recycle <- function(v, n) {
  if (length(v) == n) return(v)
  if (length(v) == 1L) return(rep(v, n))
  if (length(v) > n) return(v[seq_len(n)])
  c(v, rep(v[length(v)], n - length(v)))
}

# Intensity-weighted center of the volume, in RAS. Negative intensities (e.g.
# CT air around -1000) are clamped to 0 so the center reflects tissue/bone.
# Falls back to the geometric (field-of-view) center if the volume is empty.
center_of_mass <- function(vol, vox2ras) {
  d <- dim(vol)
  m <- as.double(vol)
  m[!is.finite(m) | m < 0] <- 0
  dim(m) <- d
  xm <- .rowSums(m, d[1], d[2] * d[3])           # mass per x slice
  yz <- .colSums(m, d[1], d[2] * d[3])           # mass per (y, z), summed over x
  dim(yz) <- c(d[2], d[3])
  ym <- rowSums(yz)
  zm <- colSums(yz)
  tot <- sum(xm)
  if (!is.finite(tot) || tot <= 0) {
    ijk <- (d - 1) / 2                            # geometric center fallback
  } else {
    ijk <- c(
      sum((0:(d[1] - 1)) * xm),
      sum((0:(d[2] - 1)) * ym),
      sum((0:(d[3] - 1)) * zm)) / tot
  }
  as.numeric(vox2ras %*% c(ijk, 1))[1:3]
}

# Deformable (symmetric diffeomorphic) stage; see src/reg_syn.cpp.
# fixed_list / moving_list are lists of co-registered channels (one or more);
# metrics has one entry per channel. weights are normalized contributions.
register_syn <- function(fixed_list, fixed_dim, fixed_vox2ras,
                          moving_list, moving_dim, moving_vox2ras,
                          metrics, weights, init_transform,
                          shrink_factors, smoothing_sigmas, iterations,
                          interp_codes, field_sigma, seed,
                          fixed_mask = NULL, moving_mask = NULL,
                          fixed_points = NULL, moving_points = NULL, points_weight = 0.5,
                          grad_step = 0.2, total_sigma = 0, cc_radius = 2L,
                          verbose = FALSE) {
  # Mattes drives the linear stage; the deformable refinement uses local CC
  # (ANTs SyN's metric) unless mean squares was requested, applied per channel.
  deform_metric <- ifelse(metrics == "meansquares", "meansquares", "cc")
  fixed_dbl <- lapply(fixed_list, as.double)
  moving_dbl <- lapply(moving_list, as.double)
  # NULL point sets become 0 x 3 matrices (the C++ side treats N = 0 as "none").
  fixed_points_mat <- if (is.null(fixed_points)) matrix(numeric(0), 0L, 3L) else as.matrix(fixed_points)
  moving_points_mat <- if (is.null(moving_points)) matrix(numeric(0), 0L, 3L) else as.matrix(moving_points)
  register_syn_cpp(
    fixedList = fixed_dbl, fixedDim = fixed_dim, fixedVox2Ras = fixed_vox2ras,
    movingList = moving_dbl, movingDim = moving_dim, movingVox2Ras = moving_vox2ras,
    affineTransform = init_transform, metricNames = deform_metric,
    weights = as.double(weights),
    shrinkFactors = shrink_factors, smoothingSigmas = smoothing_sigmas,
    iterations = iterations, gradStep = as.double(grad_step),
    flowSigma = as.double(field_sigma),   # one entry per level (recycled upstream)
    totalSigma = as.double(total_sigma),
    ccRadius = as.integer(cc_radius),
    interpCodes = as.integer(interp_codes),
    fixedMask = fixed_mask, movingMask = moving_mask,
    fixedPoints = fixed_points_mat, movingPoints = moving_points_mat,
    pointsWeight = as.double(points_weight),
    verbose = isTRUE(verbose))
}
