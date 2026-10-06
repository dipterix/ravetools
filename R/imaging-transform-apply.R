#' @title Apply a registration (linear and/or deformable) to a volume or to points
#' @description
#' \code{apply_transform3d_volume} resamples a 3D volume (or a 4D stack of
#' frames, each frame treated identically) through a chain of transforms onto a
#' reference grid; \code{apply_transform3d_points} maps \verb{RAS} coordinates
#' through the same chain. The chain is either a registration result (from
#' \code{\link{register_volume3d}} or \code{\link{load_registration}}), whose
#' mapping direction is chosen with \code{direction}, or an \pkg{'ANTs'}-style
#' list of transforms applied in the given order. Both functions stream every
#' voxel or point through the chain in \verb{C++}; no composite deformation
#' field is ever materialized. These functions supersede
#' \code{\link{apply_transform3d}}, which only handles a single linear transform.
#'
#' @section Registration objects:
#' A \code{ravetools_register_volume3d} object holds the linear transform
#' \eqn{A} (target \verb{RAS} to source \verb{RAS}) and, for \code{"syn"}
#' registrations, the forward and inverse displacement fields \eqn{u_f} and
#' \eqn{u_i}, both defined on the target grid (\verb{RAS} millimeters). The
#' warped image returned by \code{\link{register_volume3d}} samples the source
#' at \eqn{A(r + u_f(r))} for every target voxel \eqn{r}, and \code{direction}
#' selects which of the following maps is evaluated:
#'
#' \tabular{llll}{
#'   \strong{call} \tab \strong{maps} \tab \strong{composite} \tab \strong{default output grid} \cr
#'   volume, \code{"forward"} \tab a source-space image onto the target grid \tab \eqn{A(r + u_f(r))} \tab target grid \cr
#'   volume, \code{"inverse"} \tab a target-space image onto the source grid \tab \eqn{q + u_i(q)}, \eqn{q = A^{-1} r} \tab source grid \cr
#'   points, \code{"forward"} \tab source \verb{RAS} to target \verb{RAS} \tab \eqn{q + u_i(q)}, \eqn{q = A^{-1} p} \tab (none) \cr
#'   points, \code{"inverse"} \tab target \verb{RAS} to source \verb{RAS} \tab \eqn{A(p + u_f(p))} \tab (none) \cr
#' }
#'
#' Rigid and \verb{affine} registrations behave the same way with
#' \eqn{u \equiv 0}. With the same \code{direction}, a volume and a set of
#' points move the same way (\code{"forward"}: from source to target space),
#' but the composites differ: an image is pulled onto its new grid by looking
#' up where each output voxel comes from, whereas a point is pushed to where it
#' lands, so the forward point map is the inverse of the forward volume
#' composite. The inverse field is accurate wherever the forward map is
#' one-to-one inside the target grid; at the border of an unmasked
#' deformation, where the field may fold or push points out of the grid,
#' \eqn{u_i} has no exact solution and points mapped there are approximate. The
#' \code{vox2ras} of \code{volume} defaults to the argument, then to the
#' array's \code{"vox2ras"} attribute, then to the registration geometry (the
#' source grid for \code{"forward"}, the target grid for \code{"inverse"});
#' \code{reference_dim} and \code{reference_vox2ras} default to the geometry's
#' target grid for \code{"forward"} and source grid for \code{"inverse"}. A
#' manifest written before the source dimension was recorded
#' (\code{SourceDim}) still loads, but then \code{reference_dim} must be given
#' for the \code{"inverse"} direction. \code{invert = TRUE} is an error for
#' registration objects: choose the mapping with \code{direction}.
#'
#' @section Transform lists:
#' A list in \code{antsApplyTransforms -t} order. Each element is one of a
#' \eqn{4\times 4} (or \eqn{3\times 4}) \verb{RAS} matrix mapping fixed to
#' moving coordinates (such as \code{$transform} or
#' \code{\link{read_ants_transform}}), a \code{(nx, ny, nz, 3)} displacement
#' array in \verb{RAS} millimeters carrying a \code{"vox2ras"} attribute (such
#' as \code{$forward_field} or \code{\link{read_ants_warp}}), or a file path
#' (\code{.mat} read with \code{\link{read_ants_transform}}, \code{.nii} /
#' \code{.nii.gz} read with \code{\link{read_ants_warp}}). The composite is
#' \eqn{C = T_n \circ \cdots \circ T_1}: the \emph{first} element is applied
#' first to a point, exactly as \pkg{'ANTs'} does, so transform lists stored by
#' \pkg{'ANTs'}-based pipelines work unchanged. Volumes are resampled as
#' \eqn{out(r) = in(C(r))} and points as \eqn{C(p)}; a registration's
#' \code{list(forward_field, transform)} therefore reproduces the
#' \code{"forward"} volume warp, and \code{list(solve(transform), inverse_field)}
#' (or \code{list(transform, inverse_field)} with
#' \code{invert = c(TRUE, FALSE)}) the \code{"inverse"} one. \code{invert} is
#' recycled to the list length and may only be \code{TRUE} for matrices; unlike
#' \code{ANTsPy}, no inversion is ever inferred from the list layout. A bare
#' matrix or field is accepted as a one-element list. \code{reference_dim} and
#' \code{reference_vox2ras} are required with a list, and \code{direction} must
#' not be given.
#'
#' @section Sampling rules:
#' Displacement fields are sampled with \verb{trilinear} interpolation in
#' their own voxel grid; within the half-voxel margin around their outer nodes
#' the edge value is used, and beyond that extent they contribute zero
#' displacement (the rule of the \pkg{'ITK'} displacement-field transform), so
#' a point leaving the field is still carried by the remaining stages.
#' Volumes are sampled with the same
#' rounding and out-of-bounds rules as \code{\link{resample_3d_volume}} and
#' \code{\link{apply_transform3d}}: \code{'nearest'} rounds the continuous
#' voxel coordinate and requires the index to fall inside the volume,
#' \code{'trilinear'} requires the coordinate to lie within the voxel-center
#' box; voxels that do not get \code{na_fill}.
#'
#' @param volume a 3D array to resample (integer and logical arrays are
#' converted to double), or a 4D array whose frames are all transformed the
#' same way; trailing dimensions of size one beyond the fourth are dropped
#' @param transforms a \code{ravetools_register_volume3d} object, or a list of
#' transforms in \pkg{'ANTs'} order (see the sections below)
#' @param vox2ras the \eqn{4\times 4} voxel-to-\verb{RAS} matrix of
#' \code{volume}; default \code{NULL} looks for the array's \code{"vox2ras"}
#' attribute, then at the registration geometry. This and every other matrix
#' involved (reference grid, field grids, transforms) must be finite and
#' non-singular; a singular one is an error
#' @param reference_dim,reference_vox2ras the output grid: its dimension
#' (length 3) and voxel-to-\verb{RAS} matrix; required for a transform list,
#' defaulting to the registration's target (\code{"forward"}) or source
#' (\code{"inverse"}) grid otherwise
#' @param direction which mapping of a registration object to evaluate (see
#' the table); only valid with a registration object
#' @param invert logical, recycled to the list length: invert the
#' corresponding matrix before use; only valid with a transform list, and only
#' for matrices
#' @param interpolation \code{'trilinear'} (default) or \code{'nearest'} (for
#' labels or masks)
#' @param na_fill value for output voxels that fall outside \code{volume};
#' default \code{0}
#' @param points an \code{N x 3} matrix, a data frame whose first three
#' columns are \code{x, y, z}, or a length-3 vector, in \verb{RAS} millimeters;
#' rows with a missing or non-finite coordinate give \code{NA} rows
#' @returns \code{apply_transform3d_volume} returns the resampled volume with
#' dimension \code{reference_dim} (plus the frame dimension for 4D input) and
#' a \code{"vox2ras"} attribute equal to \code{reference_vox2ras}.
#' \code{apply_transform3d_points} returns an \code{N x 3} double matrix.
#' @seealso \code{\link{register_volume3d}}, \code{\link{save_registration}},
#' \code{\link{read_ants_transform}}, \code{\link{read_ants_warp}},
#' \code{\link{apply_transform3d}} (superseded, single linear transform only)
#' @examples
#'
#' # a toy registration: a blob and its shifted copy
#' nd <- c(24, 24, 24)
#' vox2ras <- diag(4); vox2ras[1:3, 4] <- -12
#' blob <- function(cx, cy, cz, s = 4) {
#'   g <- expand.grid(x = 0:(nd[1]-1), y = 0:(nd[2]-1), z = 0:(nd[3]-1))
#'   array(exp(-((g$x-cx)^2 + (g$y-cy)^2 + (g$z-cz)^2) / (2*s^2)), nd)
#' }
#' target <- blob(12, 12, 12)
#' source <- blob(14, 11, 12.5)
#' reg <- register_volume3d(source, target, vox2ras, vox2ras,
#'                          type = "rigid", metric = "cc", verbose = FALSE)
#'
#' # warp the source onto the target grid (same as reg$image) ...
#' warped <- apply_transform3d_volume(source, reg, direction = "forward")
#' max(abs(warped - reg$image))
#'
#' # ... and a target-space label map back onto the source grid
#' label <- target > 0.5
#' label_src <- apply_transform3d_volume(label, reg, direction = "inverse",
#'                                       interpolation = "nearest")
#'
#' # points: "forward" maps source RAS to target RAS (like the forward
#' # volume warp, which brings source-space data into target space)
#' p_src <- rbind(c(2, -1, 0.5), c(0, 0, 0))
#' p_tgt <- apply_transform3d_points(p_src, reg, direction = "forward")
#' apply_transform3d_points(p_tgt, reg, direction = "inverse")   # round trip
#'
#' # the same linear transform as an ANTs-style list (reference grid required)
#' same <- apply_transform3d_volume(source, list(reg$transform), vox2ras = vox2ras,
#'                                  reference_dim = nd, reference_vox2ras = vox2ras)
#' max(abs(same - warped))
#'
#' @export
apply_transform3d_volume <- function(
    volume, transforms, vox2ras = NULL,
    reference_dim = NULL, reference_vox2ras = NULL,
    direction = c("forward", "inverse"), invert = FALSE,
    interpolation = c("trilinear", "nearest"), na_fill = 0) {

  caller <- "apply_transform3d_volume"
  interpolation <- match.arg(interpolation)
  direction_given <- !missing(direction)
  direction <- match.arg(direction)
  chain <- build_transform_chain(transforms, direction = direction,
                                 direction_given = direction_given,
                                 invert = invert, mode = "volume", caller = caller)

  vol <- as_volume_frames(volume, caller)

  # moving grid: argument, then attribute, then registration geometry
  if (is.null(vox2ras)) vox2ras <- attr(volume, "vox2ras")
  if (is.null(vox2ras)) vox2ras <- chain$default_vox2ras
  if (is.null(vox2ras)) {
    stop(sprintf(
      "`%s`: `vox2ras` is missing. Pass the volume's voxel-to-RAS matrix (or set a 'vox2ras' attribute on `volume`).",
      caller))
  }
  vox2ras <- as_vox2ras_matrix(vox2ras, "vox2ras", caller)

  # reference grid: arguments, then registration geometry
  if (is.null(reference_dim)) reference_dim <- chain$default_reference_dim
  if (is.null(reference_vox2ras)) reference_vox2ras <- chain$default_reference_vox2ras
  if (is.null(reference_dim) || is.null(reference_vox2ras)) {
    stop(chain$missing_reference_message)
  }
  reference_dim <- as.integer(reference_dim)
  if (length(reference_dim) < 3L || anyNA(reference_dim) || any(reference_dim[1:3] < 1L)) {
    stop(sprintf("`%s`: `reference_dim` must be three positive integers.", caller))
  }
  reference_dim <- reference_dim[1:3]
  reference_vox2ras <- as_vox2ras_matrix(reference_vox2ras, "reference_vox2ras", caller)

  interp_code <- if (interpolation == "nearest") 0L else 1L
  na_fill <- as.double(na_fill)[[1L]]

  if (chain$all_affine) {
    # Pure linear chain: fold the matrices and use the package resampler, which
    # is what apply_transform3d does (bit-identical results).
    frames <- lapply(seq_len(vol$nframes), function(f) {
      frame <- if (vol$had4d) vol$data[, , , f] else vol$data
      dim(frame) <- vol$dim3
      resample_volume_affine(frame, vox2ras, chain$folded, reference_dim,
                             reference_vox2ras, interp_code, na_fill)
    })
    out <- if (vol$had4d) {
      array(unlist(frames, use.names = FALSE), c(reference_dim, vol$nframes))
    } else {
      frames[[1L]]
    }
  } else {
    out <- apply_transform3d_volume_cpp(
      volume = vol$data, volumeDim = vol$dim_all, volumeVox2Ras = vox2ras,
      referenceDim = reference_dim, referenceVox2Ras = reference_vox2ras,
      chain = chain$stages, interpolation = interp_code, naFill = na_fill)
  }
  dim(out) <- c(reference_dim, if (vol$had4d) vol$nframes)
  attr(out, "vox2ras") <- reference_vox2ras
  out
}

#' @rdname apply_transform3d_volume
#' @export
apply_transform3d_points <- function(
    points, transforms, direction = c("forward", "inverse"), invert = FALSE) {

  caller <- "apply_transform3d_points"
  direction_given <- !missing(direction)
  direction <- match.arg(direction)
  chain <- build_transform_chain(transforms, direction = direction,
                                 direction_given = direction_given,
                                 invert = invert, mode = "points", caller = caller)

  points <- as_points_matrix(points, caller)
  out <- matrix(NA_real_, nrow(points), 3L)
  # rows with a missing or non-finite coordinate give NA rows on every path
  ok <- rowSums(!is.finite(points)) == 0L
  if (any(ok)) {
    if (chain$all_affine) {
      m <- chain$folded
      p <- points[ok, , drop = FALSE]
      out[ok, ] <- p %*% t(m[1:3, 1:3]) + rep(m[1:3, 4], each = nrow(p))
    } else {
      out[ok, ] <- apply_transform3d_points_cpp(points[ok, , drop = FALSE], chain$stages)
    }
  }
  out
}


# ---- internal helpers -------------------------------------------------------

# Coerce a 4x4 (or 3x4) voxel-to-RAS / RAS-to-RAS matrix. Every matrix that
# reaches the appliers (volume, reference and field grids, transforms) must be
# finite with an invertible 3x3 part: a singular one would silently turn a field
# into a no-op or a volume into constant output.
as_vox2ras_matrix <- function(v, name, caller) {
  v <- as.matrix(v)
  if (!is.numeric(v) || any(!is.finite(v))) {
    stop(sprintf("`%s`: `%s` must be a finite numeric 4x4 (or 3x4) matrix.", caller, name))
  }
  if (nrow(v) == 3L && ncol(v) == 4L) v <- rbind(v, c(0, 0, 0, 1))
  if (!all(dim(v) == c(4L, 4L))) {
    stop(sprintf("`%s`: `%s` must be a 4x4 (or 3x4) matrix.", caller, name))
  }
  storage.mode(v) <- "double"
  dimnames(v) <- NULL
  if (!(rcond(v[1:3, 1:3]) > 1e-10)) {
    stop(sprintf("`%s`: `%s` must be invertible (its 3x3 part is singular or nearly so).",
                 caller, name))
  }
  v
}

# A 3D volume (or 4D stack) as double, with its dimensions split out.
as_volume_frames <- function(volume, caller) {
  d <- dim(volume)
  if (is.null(d) || length(d) < 3L) {
    stop(sprintf("`%s`: `volume` must be a 3D array (or a 4D array of frames).", caller))
  }
  if (length(d) > 4L) {
    # trailing singleton dimensions beyond the fourth are dropped (as in NIfTI
    # stacks); a fourth dimension of 1 then leaves a plain 3D volume
    if (!all(d[-(1:4)] == 1L)) {
      stop(sprintf(
        "`%s`: `volume` must be a 3D array or a 4D array of frames (dimensions beyond the fourth must be 1).",
        caller))
    }
    d <- if (d[4] == 1L) d[1:3] else d[1:4]
  }
  if (any(d < 1L)) {
    stop(sprintf("`%s`: `volume` has an empty dimension.", caller))
  }
  if (!is.numeric(volume) && !is.logical(volume)) {
    stop(sprintf("`%s`: `volume` must be numeric or logical.", caller))
  }
  storage.mode(volume) <- "double"
  dim(volume) <- d
  had4d <- length(d) == 4L
  list(data = volume, dim3 = d[1:3], dim_all = as.integer(d),
       nframes = if (had4d) d[4] else 1L, had4d = had4d)
}

# Points as an N x 3 double matrix without dimnames.
as_points_matrix <- function(points, caller) {
  if (is.data.frame(points)) {
    if (ncol(points) < 3L) {
      stop(sprintf("`%s`: `points` must have 3 columns (x, y, z in RAS).", caller))
    }
    points <- as.matrix(points[, 1:3])
  } else if (is.null(dim(points))) {
    if (length(points) != 3L) {
      stop(sprintf("`%s`: `points` must be an N x 3 matrix, a data frame, or a length-3 vector.", caller))
    }
    points <- matrix(points, 1L, 3L)
  }
  points <- as.matrix(points)
  if (ncol(points) != 3L) {
    stop(sprintf("`%s`: `points` must have 3 columns (x, y, z in RAS).", caller))
  }
  if (!is.numeric(points) && !is.logical(points)) {
    stop(sprintf("`%s`: `points` must be numeric.", caller))
  }
  storage.mode(points) <- "double"
  dimnames(points) <- NULL
  points
}

# Resample a 3D double volume through a single RAS affine (fixed -> moving)
# with the package resampler; this is exactly what apply_transform3d does.
# interp_code: 0 nearest, 1 trilinear, 2 cubic B-spline.
resample_volume_affine <- function(volume, vox2ras, transform, reference_dim,
                                   reference_vox2ras, interp_code, na_fill) {
  new_vox_to_world <- transform %*% reference_vox2ras
  storage.mode(na_fill) <- "double"
  re <- resample3D(
    arrayDim = as.integer(reference_dim[1:3]),
    fromArray = volume,
    newVoxToWorldTransposed = t(new_vox_to_world),
    oldVoxToWorldTransposed = t(vox2ras),
    na = na_fill,
    interpolation = as.integer(interp_code))
  re[[1]]
}

# Stage constructors (the list layout read by the C++ side).
affine_stage <- function(m) {
  list(type = "affine", matrix = m)
}

field_stage <- function(field, vox2ras, caller, what) {
  d <- dim(field)
  if (length(d) != 4L || d[4] != 3L) {
    stop(sprintf("`%s`: %s must be a (nx, ny, nz, 3) displacement array.", caller, what))
  }
  if (!is.numeric(field)) {
    stop(sprintf("`%s`: %s must be numeric.", caller, what))
  }
  if (!is.double(field)) storage.mode(field) <- "double"
  list(type = "field", field = field, dim = as.integer(d[1:3]),
       vox2ras = as_vox2ras_matrix(vox2ras, sprintf("vox2ras of %s", what), caller))
}

# Read one transform file as a matrix or a field.
read_transform_file <- function(path, caller) {
  lp <- tolower(path)
  if (grepl("\\.mat$", lp)) return(read_ants_transform(path))
  if (grepl("\\.nii(\\.gz)?$", lp)) return(read_ants_warp(path))
  stop(sprintf("`%s`: unsupported transform file '%s' (expected .mat, .nii or .nii.gz).",
               caller, path))
}

# Convert the i-th element of a transform list into a stage.
stage_from_element <- function(el, i, inv, caller) {
  label <- sprintf("element %d", i)
  if (is.character(el)) {
    if (length(el) != 1L || is.na(el)) {
      stop(sprintf("`%s`: %s must be a single file path.", caller, label))
    }
    label <- sprintf("element %d ('%s')", i, el)
    el <- read_transform_file(el, caller)
  }
  if (inherits(el, "ravetools_register_volume3d")) {
    stop(sprintf(
      "`%s`: %s is a registration object; a transform list may only contain matrices, fields or file paths (use `list(reg$forward_field, reg$transform)` to chain its parts).",
      caller, label))
  }
  d <- dim(el)
  if (is.numeric(el) && length(d) == 2L &&
      ((d[1] == 4L && d[2] == 4L) || (d[1] == 3L && d[2] == 4L))) {
    m <- as_vox2ras_matrix(el, label, caller)
    if (inv) m <- solve(m)
    return(affine_stage(m))
  }
  if (is.numeric(el) && length(d) == 4L && d[4] == 3L) {
    if (inv) {
      stop(sprintf(
        "`%s`: `invert` can only be TRUE for affine matrices (%s is a displacement field).",
        caller, label))
    }
    v <- attr(el, "vox2ras")
    if (is.null(v)) {
      stop(sprintf("`%s`: %s is a displacement field without a 'vox2ras' attribute.", caller, label))
    }
    return(field_stage(el, v, caller, label))
  }
  stop(sprintf(
    "`%s`: %s is not a 4x4 (or 3x4) matrix, a (nx, ny, nz, 3) displacement array, or a file path.",
    caller, label))
}

# Fold consecutive affine stages; returns the chain description used by both
# appliers.
finish_chain <- function(stages, default_vox2ras = NULL, default_reference_dim = NULL,
                         default_reference_vox2ras = NULL, missing_reference_message = "") {
  all_affine <- all(vapply(stages, function(s) identical(s$type, "affine"), logical(1)))
  folded <- NULL
  if (all_affine) {
    folded <- diag(4)
    for (s in stages) folded <- s$matrix %*% folded     # T_n %*% ... %*% T_1
    if (length(stages) == 1L) folded <- stages[[1L]]$matrix
  }
  list(stages = stages, all_affine = all_affine, folded = folded,
       default_vox2ras = default_vox2ras,
       default_reference_dim = default_reference_dim,
       default_reference_vox2ras = default_reference_vox2ras,
       missing_reference_message = missing_reference_message)
}

# Stage chain for a registration object, following the direction table in the
# documentation: volume/forward and points/inverse evaluate A(p + u_f(p));
# volume/inverse and points/forward evaluate q + u_i(q), q = A^-1 p.
registration_chain <- function(reg, direction, mode, caller) {
  A <- reg$transform
  if (is.null(A)) {
    stop(sprintf("`%s`: the registration object has no `transform` (linear part).", caller))
  }
  A <- as_vox2ras_matrix(A, "transform", caller)
  geom <- reg$geometry %||% list()
  fwd <- reg$forward_field
  inv <- reg$inverse_field
  deformable <- !is.null(fwd) || !is.null(inv) ||
    isTRUE(reg$type %in% c("syn", "syn_only"))

  one_field <- function(f, what) {
    if (is.null(f)) {
      stop(sprintf(
        "`%s`: the registration is deformable but has no `%s`, which direction = '%s' needs.",
        caller, what, direction))
    }
    v <- attr(f, "vox2ras") %||% geom$target_vox2ras
    if (is.null(v)) {
      stop(sprintf("`%s`: `%s` has no 'vox2ras' attribute and the geometry has no target grid.",
                   caller, what))
    }
    field_stage(f, v, caller, sprintf("`%s`", what))
  }

  forward_map <- (mode == "volume") == (direction == "forward")
  if (forward_map) {
    stages <- c(if (deformable) list(one_field(fwd, "forward_field")), list(affine_stage(A)))
  } else {
    stages <- c(list(affine_stage(solve(A))), if (deformable) list(one_field(inv, "inverse_field")))
  }

  if (direction == "forward") {
    finish_chain(stages,
                 default_vox2ras = geom$source_vox2ras,
                 default_reference_dim = geom$target_dim,
                 default_reference_vox2ras = geom$target_vox2ras,
                 missing_reference_message = sprintf(
                   "`%s`: the registration geometry has no target grid; supply `reference_dim` and `reference_vox2ras`.",
                   caller))
  } else {
    msg <- if (is.null(geom$source_dim)) {
      sprintf(
        "`%s`: the registration geometry has no `source_dim` (manifests written before 'SourceDim' was recorded lack it), so `reference_dim` must be supplied for direction = 'inverse'.",
        caller)
    } else {
      sprintf("`%s`: the registration geometry has no source grid; supply `reference_dim` and `reference_vox2ras`.",
              caller)
    }
    finish_chain(stages,
                 default_vox2ras = geom$target_vox2ras,
                 default_reference_dim = geom$source_dim,
                 default_reference_vox2ras = geom$source_vox2ras,
                 missing_reference_message = msg)
  }
}

# Normalize `transforms` (registration object, list, bare matrix / field /
# paths) into a stage chain plus the grid defaults.
build_transform_chain <- function(transforms, direction, direction_given, invert, mode, caller) {
  if (inherits(transforms, "ravetools_register_volume3d")) {
    if (!isFALSE(invert)) {
      stop(sprintf(
        "`%s`: `invert` must be FALSE for a registration object; choose the mapping with `direction` instead.",
        caller))
    }
    return(registration_chain(transforms, direction, mode, caller))
  }
  if (direction_given) {
    stop(sprintf(
      "`%s`: `direction` only applies to a registration object; a transform list is applied in the given order (use `invert` to invert matrices).",
      caller))
  }
  if (is.character(transforms)) {
    transforms <- as.list(transforms)
  } else if (!is.list(transforms)) {
    transforms <- list(transforms)
  }
  n <- length(transforms)
  if (n == 0L) {
    stop(sprintf("`%s`: `transforms` is empty.", caller))
  }
  invert <- as.logical(invert)
  if (length(invert) == 1L) invert <- rep(invert, n)
  if (length(invert) != n || anyNA(invert)) {
    stop(sprintf(
      "`%s`: `invert` must be a logical vector of length 1 or length(transforms) (%d), without NA.",
      caller, n))
  }
  stages <- lapply(seq_len(n), function(i) {
    stage_from_element(transforms[[i]], i, invert[i], caller)
  })
  finish_chain(stages, missing_reference_message = sprintf(
    "`%s`: `reference_dim` and `reference_vox2ras` are required when `transforms` is a list of transforms.",
    caller))
}

# Numerically accurate inverse of a dense displacement field by damped
# fixed-point iteration (see src/reg_apply.h); used by register_volume3d for
# its `inverse_field`. `field` is (nx, ny, nz, 3) in RAS millimeters on the grid
# described by `vox2ras`; tolerances are in voxel units. Returns the inverse on
# the same grid with a "vox2ras" attribute and an "inversion" attribute holding
# the iteration count and the final residuals.
invert_displacement_field <- function(field, vox2ras = NULL, max_iterations = 50L,
                                      mean_tolerance = 1e-5, max_tolerance = 1e-3) {
  caller <- "invert_displacement_field"
  if (is.null(vox2ras)) vox2ras <- attr(field, "vox2ras")
  if (is.null(vox2ras)) {
    stop(sprintf("`%s`: `vox2ras` is required (argument or 'vox2ras' attribute).", caller))
  }
  vox2ras <- as_vox2ras_matrix(vox2ras, "vox2ras", caller)
  d <- dim(field)
  if (length(d) != 4L || d[4] != 3L) {
    stop(sprintf("`%s`: `field` must be a (nx, ny, nz, 3) array.", caller))
  }
  if (!is.double(field)) storage.mode(field) <- "double"
  res <- invert_displacement_field_cpp(
    field = field, dim = as.integer(d[1:3]), vox2ras = vox2ras,
    maxIterations = as.integer(max_iterations)[[1L]],
    meanTolerance = as.double(mean_tolerance)[[1L]],
    maxTolerance = as.double(max_tolerance)[[1L]])
  out <- res$field
  dim(out) <- d
  attr(out, "vox2ras") <- vox2ras
  attr(out, "inversion") <- list(iterations = res$iterations,
                                 mean_residual = res$mean_residual,
                                 max_residual = res$max_residual)
  out
}
