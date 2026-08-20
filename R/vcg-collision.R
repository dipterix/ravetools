# Resolve one side of `vcg_detect_collision` into the layout the C++ core
# expects: a 3 x n coordinate matrix plus, in mesh mode, a 0-based 3 x m face
# matrix. Returns the resolved mode alongside, since "auto" is decided here.
resolve_collision_input <- function(obj, mode, arg) {
  is_matrix_like <- is.matrix(obj) || (is.array(obj) && length(dim(obj)) == 2L)

  if (mode == "auto" && !is_matrix_like) {
    mesh <- ensure_mesh3d(obj)
    if (length(mesh$it)) {
      mode <- "mesh"
    }
  }

  if (mode == "auto") {
    # Check points or segments
    mode <- "points"
  }

  if (mode == "mesh") {
    if (is_matrix_like) {
      stop(sprintf(
        "`vcg_detect_collision`: `%s` is a matrix, which carries no faces; `mode_%s = \"mesh\"` needs a `mesh3d` object. Use \"points\" or \"segments\" instead.",
        arg, arg))
    }
    surface <- meshintegrity(mesh = ensure_mesh3d(obj), facecheck = TRUE)
    return(list(
      mode = mode,
      vb = surface$vb[1:3, , drop = FALSE],
      it = surface$it - 1L
    ))
  }

  vb <- as_point_cloud_matrix(obj, "vcg_detect_collision")

  # Normalize partial NA / non-finite rows to whole NA columns so the C++ side
  # only ever sees a column that is either wholly usable or a separator.
  incomplete <- !is.finite(vb[1L, ]) | !is.finite(vb[2L, ]) | !is.finite(vb[3L, ])
  if (any(incomplete)) {
    vb[, incomplete] <- NA_real_
  }

  list(mode = mode, vb = vb, it = NULL)
}


#' @title Detect collisions between two geometries
#' @description
#' Reports, for every element of \code{y}, whether it comes within
#' \code{radius} of \code{x}, together with the exact minimum distance. Each
#' side may be a point cloud, a set of connected line segments (for example
#' diffusion streamlines or an electrode shaft), or a triangular mesh, so the
#' function covers all nine pairings with one call.
#'
#' The test is exact: it produces neither false negatives nor false positives,
#' and \code{distance} is the true minimum distance rather than a
#' vertex-sampled approximation.
#'
#' @param x the geometry that gets indexed; a matrix with \code{n} rows and 2
#' or 3 columns, or a surface (see \code{\link{ensure_mesh3d}})
#' @param y the geometry that gets queried; results are aligned to this side.
#' Same accepted types as \code{x}
#' @param mode_x,mode_y how to interpret \code{x} and \code{y}: \code{'points'}
#' treats each row as an isolated point, \code{'segments'} chains consecutive
#' rows together (\code{p1-p2}, \code{p2-p3}, ...), and \code{'mesh'} uses
#' the triangular faces. The default \code{'auto'} resolves to \code{'mesh'}
#' for a surface that has faces and \code{'points'} otherwise, so a matrix of
#' connected line segments must ask for \code{'segments'} explicitly
#' @param radius distance tolerance; a collision is reported when the minimum
#' distance is at most \code{radius}. Defaults to \code{0}, meaning literal
#' contact; a positive tolerance is usually what is wanted
#' @param early_stop whether to stop as soon as a collision is found, scanning
#' independently within each group of \code{y} delimited by \code{NA} rows. A
#' matrix carrying no \code{NA} rows is a single group, so a point cloud stops
#' after one hit overall rather than per row; likewise \code{mode_y = 'mesh'}
#' stops at the first colliding face. Skipped elements report \code{hit = NA},
#' meaning \emph{not evaluated} rather than \emph{no collision}, and the
#' collision reported is the first in row order, not the closest. Use
#' \code{FALSE} (the default) whenever the per-element profile or a minimum
#' distance is wanted
#' @param include_interior whether geometry lying strictly inside a closed
#' \code{x} surface counts as a collision even when it never comes within
#' \code{radius} of the surface itself; requires \code{mode_x = 'mesh'}.
#' Default is \code{FALSE}
#'
#' @returns A list of three vectors, each aligned one-to-one with the rows of
#' \code{y} (or with the faces of \code{y} when \code{mode_y = 'mesh'}):
#' \describe{
#'   \item{\code{hit}}{logical. \code{NA} marks an element with nothing to
#'     test: a separator row, or the final vertex of each group when
#'     \code{mode_y = 'segments'}. When \code{early_stop = TRUE}, \code{NA}
#'     also marks elements the scan never reached}
#'   \item{\code{distance}}{the exact minimum distance where \code{hit} is
#'     \code{TRUE}, and \code{NA} everywhere else. Distances beyond
#'     \code{radius} are never computed, so a \code{FALSE} element carries no
#'     distance}
#'   \item{\code{index}}{1-based index of the closest element of \code{x}
#'     within \code{radius}: a row when \code{mode_x} is \code{'points'}, the
#'     row where the segment starts when \code{'segments'}, or a face column
#'     when \code{'mesh'}. \code{NA} where there is no collision}
#' }
#'
#' @details
#' \code{x} is indexed once into a uniform spatial grid and \code{y} is
#' streamed against it, so put the larger or repeatedly-reused geometry in
#' \code{x} and the geometry you want per-element answers about in \code{y}.
#'
#' Several chains share one matrix and are delimited by rows of \code{NA}; any
#' row that is not fully finite is treated as a separator. This lets a whole
#' bundle of streamlines be passed as a single matrix:
#'
#' \preformatted{
#'   rbind(
#'     c(0, 0, 0), c(1, 0, 0), c(2, 0, 0),   # first streamline
#'     c(NA, NA, NA),                        # separator
#'     c(0, 5, 0), c(1, 5, 0)                # second streamline
#'   )
#' }
#'
#' Those separators also define the unit that \code{early_stop} works on: it
#' stops once per group, so a bundle of streamlines yields at most one hit per
#' streamline.
#'
#' When \code{mode_y = 'segments'}, row \code{i} of the result describes the
#' segment running from row \code{i} to row \code{i + 1}.
#'
#' Both geometries must already share one coordinate space. Surface
#' coordinates, volume \verb{IJK} indices, and scanner coordinates are all
#' different spaces, and no transform is applied here.
#'
#' \code{include_interior} relies on ray casting and therefore needs \code{x}
#' to be watertight with coherently oriented faces; see
#' \code{\link{vcg_fix_defects}} to repair a surface that is not.
#'
#' Multi-threading follows \code{\link{ravetools_threads}}. The interior test
#' always runs single-threaded.
#'
#' @inheritSection ensure_mesh3d Coercing Surface Inputs
#'
#' @examples
#'
#' library(ravetools)
#'
#' # A spherical region of interest
#' roi <- vcg_sphere()
#'
#' # Two streamlines in one matrix, separated by an NA row: the first passes
#' # through the sphere, the second stays well clear of it
#' streamlines <- rbind(
#'   cbind(seq(-3, 3, by = 0.5), 0, 0),
#'   c(NA, NA, NA),
#'   cbind(seq(-3, 3, by = 0.5), 5, 0)
#' )
#'
#' result <- vcg_detect_collision(
#'   x = roi,
#'   y = streamlines,
#'   mode_y = "segments",
#'   radius = 0.1
#' )
#'
#' # Which segments touch the sphere, and how close do they get?
#' data.frame(
#'   hit = result$hit,
#'   distance = round(result$distance, 4)
#' )
#'
#' # A point cloud against the same region, asking only for proximity
#' points <- rbind(c(0, 0, 0), c(1.05, 0, 0), c(10, 10, 10))
#' vcg_detect_collision(roi, points, radius = 0.1)$hit
#'
#' # The centre of the sphere is far from its surface, so it only counts as a
#' # collision when the interior is included
#' vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1)$hit
#' vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1,
#'                      include_interior = TRUE)$hit
#'
#' @export
vcg_detect_collision <- function(
    x, y,
    mode_x = c("auto", "points", "segments", "mesh"),
    mode_y = c("auto", "points", "segments", "mesh"),
    radius = 0,
    early_stop = FALSE,
    include_interior = FALSE) {

  mode_x <- match.arg(mode_x)
  mode_y <- match.arg(mode_y)

  radius <- as.double(radius)
  stopifnot2(
    length(radius) == 1L && !is.na(radius) && is.finite(radius) && radius >= 0,
    msg = "`vcg_detect_collision`: `radius` must be a single finite non-negative number."
  )

  early_stop <- isTRUE(as.logical(early_stop))
  include_interior <- isTRUE(as.logical(include_interior))

  side_x <- resolve_collision_input(x, mode_x, "x")
  side_y <- resolve_collision_input(y, mode_y, "y")

  if (include_interior && side_x$mode != "mesh") {
    stop("`vcg_detect_collision`: `include_interior = TRUE` needs `x` to be a closed surface, so `mode_x` must resolve to \"mesh\" (it resolved to \"", side_x$mode, "\").")
  }

  mode_code <- c(points = 0L, segments = 1L, mesh = 2L)

  vcgDetectCollision(
    x_ = side_x$vb,
    x_it_ = side_x$it,
    y_ = side_y$vb,
    y_it_ = side_y$it,
    mode_x = mode_code[[side_x$mode]],
    mode_y = mode_code[[side_y$mode]],
    radius = radius,
    early_stop = early_stop,
    include_interior = include_interior
  )
}
