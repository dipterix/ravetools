#' @title Outer (envelope) surface of a cortical surface
#' @description
#' Computes the smooth surface that wraps a \verb{pial} surface without
#' dipping into the \verb{sulci}: the surface \verb{FreeSurfer} writes as
#' \verb{?h.pial-outer-smoothed} when it computes the local
#' \verb{gyrification} index (\code{mris_compute_lgi}, run by
#' \code{recon-all -localGI}), re-implemented here without 'MATLAB'.
#'
#' The envelope is not a convex hull. It is a morphological closing of the
#' volume inside the \verb{pial} surface: a dilation followed by an erosion
#' with a disk or ball of radius \code{close_radius}, which fills every gap the
#' element cannot enter. The \verb{sulci}, which are narrower than the
#' element, are bridged, while concavities broader than the element are
#' followed. A convex hull is the limit of an infinitely large element.
#'
#' @details
#' With \code{method = "freesurfer"} (default), the steps reproduce
#' \code{mris_compute_lgi} as it is coded:
#' \enumerate{
#'   \item Fill the lattice points that lie inside \code{surface}, as
#'     \code{mris_fill -c -r 1} does. The lattice is 1 millimeter apart and
#'     aligned with integer surface coordinates (the conformed volume). The
#'     inside test is evaluated at the lattice points shifted by \code{delta}
#'     (less than half a millimeter per axis): \code{mris_fill} fills a grid
#'     anchored at the smallest vertex coordinate and then re-samples it onto
#'     the conformed grid with nearest-neighbor interpolation.
#'   \item Grow the filled volume by one voxel toward the right (\verb{+x})
#'     and superior (\verb{+z}) sides only. This is what the "soft Gaussian
#'     filter" of \code{make_outer_surface.m} (a 2-by-2 kernel followed by a
#'     threshold) amounts to.
#'   \item Close the volume with a disk of radius \code{close_radius}
#'     millimeters, independently in every coronal slice. The structuring
#'     element \code{strel('ball', 15, 15)} in \code{make_outer_surface.m} is
#'     two-dimensional, so nothing is bridged along the anterior-posterior
#'     axis.
#'   \item Extract the surface of the closed 0/255 volume with marching cubes
#'     at level 100, as \code{isosurface} does in \code{make_outer_surface.m}.
#'   \item Keep the largest connected component, as
#'     \code{mris_extract_main_component} does.
#'   \item Apply \code{smooth_iterations} rounds of include-self neighbor
#'     averaging, as \code{mris_smooth -nw -n 30} does (see
#'     \code{\link{mris_smooth}}).
#' }
#' On both hemispheres of the \verb{N27} template, every vertex of the
#' surfaces \verb{FreeSurfer} generated lies within 0.9 millimeters of the
#' result (mean nearest-vertex distance 0.34 millimeters).
#'
#' With \code{method = "isotropic"}, the volume is closed with a
#' three-dimensional ball of radius \code{close_radius} instead, with no
#' one-sided growth, no \code{delta} shift, and marching cubes at (just below)
#' the midpoint level. This is the procedure the \verb{FreeSurfer}
#' documentation describes. The result does not depend on how the head is
#' oriented in the volume, but it differs from the \verb{FreeSurfer} output by
#' up to a few millimeters.
#'
#' @param surface \verb{pial} surface of one hemisphere, in \verb{FreeSurfer}
#' surface (\verb{tkrRAS}) coordinates, for example read with
#' \code{freesurferformats::read.fs.surface}; see \code{\link{ensure_mesh3d}}
#' for the accepted classes. The surface must be closed, but it does not need
#' to be free of topological defects
#' @param close_radius radius of the closing element in millimeters; default
#' is \code{15}, the \code{--close_sphere_size} default of
#' \code{mris_compute_lgi}. That script documents the value as the sphere
#' diameter, but its code uses it as the radius; pass \code{7.5} for a
#' 15-millimeter diameter. A value of \code{0} skips the closing
#' @param smooth_iterations number of neighbor-averaging rounds applied to the
#' extracted surface; default is \code{30}, the \code{--smooth_iters} default
#' of \code{mris_compute_lgi}
#' @param method \code{"freesurfer"} (default) to reproduce the
#' \verb{FreeSurfer} output, or \code{"isotropic"} for a three-dimensional
#' closing; see 'Details'
#' @param verbose whether to print progress messages; default is \code{FALSE}
#'
#' @returns A triangular mesh of class \code{'mesh3d'} (with \code{vb},
#' \code{it}, and \code{normals}) in the coordinates of \code{surface}. Its
#' attribute \code{"info"} is a named list with \code{method},
#' \code{close_radius}, \code{smooth_iterations}, \code{delta} (the shift of
#' the inside test, in millimeters), \code{grid_origin} and \code{grid_dim}
#' (the corner coordinates and dimension of the working lattice),
#' \code{odd_columns} (lattice columns that crossed \code{surface} an odd
#' number of times and were left empty; non-zero for open or
#' self-intersecting input), \code{n_components} (connected components of
#' the extracted surface before the largest one is kept),
#' \code{is_closed_manifold} (see \code{\link{vcg_count_edge_defects}}), and
#' \code{euler} (the \verb{Euler} characteristic of the result; \code{2} for a
#' closed surface without handles).
#'
#' @references
#' A surface-based approach to quantify local cortical \verb{gyrification}.
#' \emph{IEEE Transactions on Medical Imaging}, 27(2), 161-170 (2008).
#'
#' @seealso \code{\link{mris_smooth}}, \code{\link{fill_surface}}
#'
#' @examples
#'
#' # Two boxes 4 millimeters apart: the envelope bridges the gap
#' volume <- array(0, c(64, 40, 40))
#' volume[10:28, 10:30, 10:30] <- 1
#' volume[33:51, 10:30, 10:30] <- 1
#' boxes <- vcg_isosurface(volume, threshold_lb = 0.5, vox_to_ras = diag(4))
#'
#' envelope <- mris_outer_surface(boxes, close_radius = 8)
#' attr(envelope, "info")$n_components
#'
#' plot_mesh_polygon(
#'   list(boxes, envelope),
#'   alpha = c(1, 0.3),
#'   col = list("gray", "red"),
#'   main = "Gray: input; red: outer surface"
#' )
#'
#'
#' @export
mris_outer_surface <- function(
    surface,
    close_radius      = 15,
    smooth_iterations = 30L,
    method            = c("freesurfer", "isotropic"),
    verbose           = FALSE
) {
  method <- match.arg(method)
  close_radius <- as.double(close_radius)
  if (length(close_radius) != 1L || !is.finite(close_radius) ||
      close_radius < 0) {
    stop("`mris_outer_surface`: `close_radius` must be a single non-negative number")
  }
  smooth_iterations <- as.integer(smooth_iterations)
  if (length(smooth_iterations) != 1L || is.na(smooth_iterations) ||
      smooth_iterations < 0L) {
    stop("`mris_outer_surface`: `smooth_iterations` must be a single non-negative integer")
  }
  verbose <- isTRUE(as.logical(verbose))
  freesurfer <- method == "freesurfer"

  mesh <- meshintegrity(surface, facecheck = TRUE)
  vertices <- mesh$vb[1:3, , drop = FALSE]

  # Working lattice: 1 millimeter apart and aligned with integer surface
  # coordinates, as the conformed volume `mris_fill -c` writes, cropped to the
  # surface plus room for the one-voxel growth, the closing, and an empty
  # border
  lo <- apply(vertices, 1L, min)
  hi <- apply(vertices, 1L, max)
  margin <- ceiling(close_radius) + 4
  grid_origin <- floor(lo) - margin
  grid_dim <- as.integer(ceiling(hi) + margin - grid_origin + 1)

  # `mris_fill` fills a grid anchored at the smallest vertex coordinate and
  # re-samples it (nearest neighbor) onto the conformed grid, so the inside
  # test ends up evaluated at the lattice points shifted by `delta`
  delta <- if (freesurfer) lo - round(lo) else c(0, 0, 0)

  if (verbose) {
    message(sprintf("Filling the surface interior on a %s lattice",
                    paste(grid_dim, collapse = " x ")))
  }
  volume <- mrisFillInteriorParity(
    vb_ = vertices - (grid_origin + delta),
    it_ = mesh$it,
    dims = grid_dim
  )
  odd_columns <- attr(volume, "odd_columns")
  if (odd_columns > 0) {
    warning(sprintf(paste(
      "`mris_outer_surface`: %d lattice columns crossed `surface` an odd",
      "number of times and were left empty; `surface` may be open or",
      "self-intersecting"), odd_columns))
  }
  if (!any(volume)) {
    stop("`mris_outer_surface`: no lattice point lies inside `surface`; is it closed?")
  }

  if (freesurfer) {
    volume <- outer_surface_grow_right_superior(volume)
  }
  if (close_radius > 0) {
    if (verbose) {
      message(sprintf("Closing with a %s of radius %g",
                      if (freesurfer) "coronal disk" else "ball", close_radius))
    }
    volume <- mrisCloseBinary(
      mask = volume,
      dims = grid_dim,
      radius = close_radius,
      slice_axis = if (freesurfer) 1L else -1L
    )
  }

  # Marching cubes on the 0/255 volume: at level 100, the vertices sit where
  # `isosurface` places them in `make_outer_surface.m` (155/255 of the way
  # from the inside voxel to the outside one). The isotropic level stays just
  # below the midpoint: a cube face with two diagonal inside corners has its
  # saddle value at exactly 127.5, where the two cubes sharing the face can
  # triangulate it inconsistently and leave cracks
  if (verbose) {
    message("Extracting the surface")
  }
  iso <- vcgIsoSurface(volume * 255L, if (freesurfer) 100 else 127)
  if (!length(iso$it)) {
    stop("`mris_outer_surface`: no surface could be extracted")
  }
  envelope <- structure(
    list(vb = iso$vb + grid_origin, it = iso$it, normals = iso$normals),
    class = c("ravetools_mesh3d", "mesh3d")
  )
  if (!checkFaceOrientation(envelope)) {
    envelope <- invertFaces(envelope)
  }

  # Keep the largest connected component; a face's vertices always belong to
  # the same component
  components <- mrisLargestComponent(
    it_ = envelope$it, nv = ncol(envelope$vb))
  keep <- components$keep
  vb <- envelope$vb[1:3, keep, drop = FALSE]
  it <- envelope$it[, keep[envelope$it[1, ]], drop = FALSE]
  it[] <- cumsum(keep)[it]

  if (verbose) {
    message(sprintf("Smoothing %d vertices for %d iterations",
                    ncol(vb), smooth_iterations))
  }
  smoothed <- mrisSmooth(
    vb_         = vb,
    it_         = it,
    niterations = smooth_iterations,
    npasses     = 1L,
    rescale     = FALSE,
    verbose     = FALSE
  )

  result <- structure(
    list(
      vb      = rbind(smoothed$vb, 1),
      it      = it,
      normals = rbind(smoothed$normals, 1)
    ),
    class = c("ravetools_mesh3d", "mesh3d")
  )
  attr(result, "info") <- list(
    method             = method,
    close_radius       = close_radius,
    smooth_iterations  = smooth_iterations,
    delta              = delta,
    grid_origin        = grid_origin,
    grid_dim           = grid_dim,
    odd_columns        = odd_columns,
    n_components       = components$n_components,
    is_closed_manifold = vcg_count_edge_defects(result)$is_closed_manifold,
    euler              = mesh_euler_characteristic(it, ncol(vb))
  )
  result
}

# `make_outer_surface.m` convolves each coronal slice of the 0/255 filled
# volume with `fspecial('gaussian', [2 2], 1)`, whose four weights are equal,
# and keeps the values above 25: a voxel is set when it or its neighbor at
# -x, -z, or (-x, -z) is set, a one-voxel growth toward +x (right) and +z
# (superior) only. `mask` is a logical array indexed by (x, y, z).
outer_surface_grow_right_superior <- function(mask) {
  d <- dim(mask)
  out <- mask
  out[-1L, , ] <- out[-1L, , ] | mask[-d[[1]], , ]
  out[, , -1L] <- out[, , -1L] | mask[, , -d[[3]]]
  out[-1L, , -1L] <- out[-1L, , -1L] | mask[-d[[1]], , -d[[3]]]
  out
}

# vertices - edges + faces of a triangle mesh (`it`: 3 x nf, one-based)
mesh_euler_characteristic <- function(it, n_vertices) {
  a <- c(it[1L, ], it[2L, ], it[3L, ])
  b <- c(it[2L, ], it[3L, ], it[1L, ])
  edge_key <- (pmin(a, b) - 1) * as.double(n_vertices) + pmax(a, b)
  n_vertices - length(unique(edge_key)) + ncol(it)
}
