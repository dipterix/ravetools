#' @title Inflate a cortical surface mesh
#' @description
#' Iteratively relaxes a closed cortical-surface mesh into a smoother, more
#' compact "inflated" shape, while tracking how far each vertex moves inward
#' along the surface normal as it goes. Returns both the inflated mesh and
#' this per-vertex depth map (\code{sulc}).
#'
#' @details
#' The implementation follows the inflation procedure described in the
#' literature (see \strong{References}):
#' \enumerate{
#'   \item Optionally rescale the mesh to a canonical surface area
#'     (110,000 \eqn{mm^2}, the normalization target used by the reference
#'     procedure).
#'   \item Outer loop: a neighborhood-averaging size \code{n_averages} is
#'     halved at each level (16, 8, ..., 0).
#'   \item Inner loop (\code{niterations} repetitions per level):
#'     \describe{
#'       \item{Distance term}{Restoring force pulling each vertex back towards
#'         its original distances to its neighbors, with coefficient
#'         \code{l_dist * sqrt(n_averages)}.}
#'       \item{Gradient averaging}{Smooth the distance-term gradient over
#'         \code{n_averages} neighborhood passes before adding the spring term.}
#'       \item{Normalized spring term}{Laplacian (neighbor-averaging) smoothing
#'         scaled by \code{sqrt(orig_area / current_area)}, with coefficient
#'         \code{l_spring_norm}; added \emph{after} gradient averaging.}
#'       \item{Momentum step}{Update vertex positions using momentum
#'         integration with a 1 mm per-step displacement cap.}
#'       \item{Depth accumulation}{Accumulate the normal-projected component
#'         of each step into the per-vertex depth map (\code{sulc}).}
#'     }
#'   \item Center the mesh, rescale it, and zero-mean the depth map (\code{sulc}).
#' }
#'
#' @param mesh triangular mesh of class \code{'mesh3d'}, typically a white-matter
#'   surface in millimeters. The mesh \strong{must} be closed, manifold, and
#'   genus-0 (watertight, no boundary or non-manifold edges, single connected
#'   component) - this is a hard precondition of the underlying algorithm. An
#'   open or defective mesh (e.g. a raw surface extracted via
#'   \code{\link{vcg_isosurface}}, which often contains small extraction
#'   cracks) will silently inflate into a distorted, "sausage"-shaped result
#'   rather than a smooth sphere-like surface, and \code{mris_inflate} will
#'   raise an error if such defects are detected.
#' @param n_averages starting number of gradient-averaging passes (outer loop);
#'   halved each level down to 0. Default \code{16}.
#' @param niterations number of inner iterations per averaging level.
#'   Default \code{10}.
#' @param l_spring_norm normalized spring term coefficient. Default \code{1.0}.
#' @param l_dist distance-preservation coefficient (before per-level scaling).
#'   Default \code{0.1}.
#' @param momentum momentum coefficient. Default \code{0.9}.
#' @param dt time step. Default \code{0.9}.
#' @param desired_rms target \verb{RMS} tangent-plane height for early stopping
#'   of each inner loop. Default \code{0.015} (mm).
#' @param scale_brain logical; whether to scale the mesh to the canonical
#'   surface area (110,000 \eqn{mm^2}) before inflation. Set \code{FALSE} only
#'   for non-standard inputs. Default \code{TRUE}.
#' @param verbose logical; print per-iteration progress. Default \code{FALSE}.
#'
#' @returns A named list:
#' \describe{
#'   \item{\code{mesh}}{Inflated surface as a \code{'mesh3d'} object with
#'     \code{vb}, \code{it}, and \code{normals}.}
#'   \item{\code{sulc}}{Numeric vector of the per-vertex depth values described
#'     above (zero-mean, in mesh units).}
#' }
#'
#' @references
#' Cortical surface-based analysis II: Inflation, flattening, and a
#' surface-based coordinate system. \emph{NeuroImage}, 9(2), 195-207 (1999).
#'
#' @examples
#' if (is_not_cran()) {
#'
#'
#' data("left_hippocampus_mask")
#' mesh <- vcg_isosurface(left_hippocampus_mask)
#'
#' # Fix defects
#' mesh <- vcg_fix_defects(mesh, verbose = TRUE)
#'
#' # Center the mesh
#' mesh$vb[1:3, ] <- mesh$vb[1:3, ] - rowMeans(mesh$vb[1:3, ])
#'
#' # Inflate the surface while keeping the node distances
#' result <- mris_inflate(mesh, n_averages = 4L, niterations = 5L,
#'                        scale_brain = FALSE, verbose = TRUE)
#'
#' # Visualize with the sulcal values
#' pal <- colorRampPalette(c("black", "gray", "red"))(128)
#' col <- pal[pmax(pmin(round(result$sulc * 10 + 64), 128), 1)]
#'
#' oldpar <- par(mfrow = c(1, 2))
#' on.exit({ par(oldpar) })
#'
#' plot(
#'   mesh, col = col,
#'   eye = c(0, 100, 0),
#'   up = c(1, 0, 0))
#'
#' plot(
#'   result$mesh, col = col,
#'   eye = c(0, 100, 0),
#'   up = c(1, 0, 0))
#'
#'
#' }
#'
#' @export
mris_inflate <- function(mesh,
                         n_averages    = 16L,
                         niterations   = 10L,
                         l_spring_norm = 1.0,
                         l_dist        = 0.1,
                         momentum      = 0.9,
                         dt            = 0.9,
                         desired_rms   = 0.015,
                         scale_brain   = TRUE,
                         verbose       = FALSE) {
  mesh <- meshintegrity(mesh, facecheck = TRUE)

  vb <- mesh$vb[1:3, , drop = FALSE]
  storage.mode(vb) <- "double"

  it <- mesh$it
  storage.mode(it) <- "integer"

  tmp <- mrisInflate(
    vb_           = vb,
    it_           = it,
    n_averages    = as.integer(n_averages)[[1L]],
    niterations   = as.integer(niterations)[[1L]],
    l_spring_norm = as.double(l_spring_norm)[[1L]],
    l_dist        = as.double(l_dist)[[1L]],
    momentum      = as.double(momentum)[[1L]],
    dt            = as.double(dt)[[1L]],
    desired_rms   = as.double(desired_rms)[[1L]],
    scale_brain   = as.logical(scale_brain)[[1L]],
    verbose       = as.logical(verbose)[[1L]]
  )

  inflated <- structure(list(
    vb      = rbind(tmp$vb, 1),
    it      = it,
    normals = rbind(tmp$normals, 1)
  ),
  class = c("ravetools_mesh3d", "mesh3d"))

  list(mesh = inflated, sulc = as.numeric(tmp$sulc))
}
