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
#' @param test_level how exhaustively \code{y} is scanned. Answers are always
#' reported per \emph{unit} of \code{y}: a point when \code{mode_y = 'points'},
#' a chain when \code{'segments'}, a face when \code{'mesh'}.
#' \code{'element'} (the default) measures every element of every unit and
#' reports the closest one; \code{'unit'} stops at the first hit inside each
#' unit and reports that one; \code{'whole'} stops at the first hit anywhere and
#' answers only whether the two geometries collide. Because a point and a face
#' each \emph{are} their own unit, \code{'element'} and \code{'unit'} differ
#' only when \code{mode_y = 'segments'}
#' @param include_interior whether geometry lying strictly inside a closed
#' \code{x} surface counts as a collision even when it never comes within
#' \code{radius} of the surface itself; requires \code{mode_x = 'mesh'}.
#' Default is \code{FALSE}
#'
#' @returns A list. \code{collide} is always present, so a plain yes-or-no
#' question can be asked at any \code{test_level}:
#' \describe{
#'   \item{\code{collide}}{a length-one logical: whether \code{x} and \code{y}
#'     collide at all. \code{FALSE} for empty input}
#'   \item{\code{hit_unit}}{logical, one entry per unit of \code{y}, omitted
#'     when \code{test_level = 'whole'}. \code{NA} marks a unit with nothing to
#'     test, and is deliberately not \code{FALSE}; see \emph{Why
#'     \code{hit_unit} reports \code{NA}} below}
#'   \item{\code{representation}}{a data frame with one row per \emph{hit}
#'     unit, so \code{nrow} equals \code{sum(hit_unit, na.rm = TRUE)} and a
#'     collision-free query gives zero rows. When
#'     \code{test_level = 'whole'} it holds at most one row}
#'   \item{\code{summary}}{what the query actually was, so a stored result can
#'     be read back without the call that made it}
#' }
#'
#' \code{summary} holds \code{radius}, \code{test_level},
#' \code{include_interior}, and one entry per side. Each side reports its
#' \code{mode} \emph{as resolved} (never \code{'auto'}), the matching
#' \code{unit_type} of \code{'point'}, \code{'line'}, or \code{'face'}, and
#' \code{n_units}: one per row for \code{'points'} with separator rows included,
#' one per chain for \code{'segments'}, one per face for \code{'mesh'}. So
#' \code{summary$y$n_units} always equals \code{length(hit_unit)}, and remains
#' the only record of the size of \code{y} when \code{test_level = 'whole'}
#' drops that vector.
#'
#' The columns of \code{representation} name the collision on both sides:
#' \describe{
#'   \item{\code{unit}}{1-based unit of \code{y}: the row itself when
#'     \code{mode_y = 'points'}, the position of the chain when
#'     \code{'segments'}, the face column when \code{'mesh'}}
#'   \item{\code{index}}{1-based element \emph{within} that unit: always
#'     \code{1} for a point, the position of the segment counted from the start
#'     of its own chain for \code{'segments'}, and the face vertex nearest
#'     \code{x} (\code{1}, \code{2}, or \code{3}) for \code{'mesh'}}
#'   \item{\code{distance}}{the exact distance from that element to \code{x}.
#'     Distances beyond \code{radius} are never computed}
#'   \item{\code{x_index}}{1-based element of \code{x} that was hit: a row when
#'     \code{mode_x} is \code{'points'}, the row where the segment starts when
#'     \code{'segments'}, or a face column when \code{'mesh'}. \code{NA} for a
#'     hit found by \code{include_interior} rather than by proximity}
#' }
#'
#' @details
#' \code{x} is indexed once into a uniform spatial grid and \code{y} is
#' streamed against it, so put the larger or repeatedly-reused geometry in
#' \code{x} and the geometry you want per-unit answers about in \code{y}.
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
#' Those separators define the units that results are reported against: with
#' \code{mode_y = 'segments'} the matrix above is two units, so
#' \code{hit_unit} has two entries whatever the \code{test_level}. With
#' \code{mode_y = 'points'} every row is its own unit, separators included, so
#' \code{hit_unit} stays aligned with the rows of \code{y} and reports
#' \code{NA} at each separator.
#'
#' Separators split the matrix into slots, and \emph{every} slot is a chain,
#' empty ones included: two adjacent separators enclose an empty chain, which
#' reports \code{NA}. The two ends are deliberately not symmetric. A separator
#' in the first row opens an empty leading chain, but a separator in the last
#' row closes the chain before it rather than opening a trailing empty one:
#'
#' \tabular{ll}{
#'   \strong{matrix} \tab \strong{chains} \cr
#'   \code{A}                \tab 1: \code{[A]} \cr
#'   \code{A NA B}           \tab 2: \code{[A] [B]} \cr
#'   \code{A NA NA B}        \tab 3: \code{[A] [] [B]} \cr
#'   \code{NA A}             \tab 2: \code{[] [A]} \cr
#'   \code{A NA}             \tab 1: \code{[A]} \cr
#' }
#'
#' That rule exists so a bundle assembled by appending a separator to every
#' chain gives exactly one unit per chain, whichever of them happen to be empty:
#'
#' \preformatted{
#'   do.call("rbind", lapply(tracts, function(line) rbind(line, NA)))
#' }
#'
#' Since \code{mode_x} plays no part in that grouping, the \code{index} of a
#' segment is always counted from the start of its own chain rather than from
#' the top of the matrix. The \code{unit} column gives the chain.
#'
#' @section Why \code{hit_unit} reports \code{NA}:
#'
#' \code{NA} means \emph{there was no geometry here to test}. \code{FALSE} means
#' \emph{geometry was here and it missed}. Three situations produce the first:
#'
#' \tabular{ll}{
#'   \strong{\code{mode_y}} \tab \strong{the unit reports \code{NA} when} \cr
#'   \code{'points'}   \tab the row is a separator \cr
#'   \code{'segments'} \tab the chain is empty (two adjacent separators) \cr
#'   \code{'segments'} \tab the chain holds one vertex, so carries no segment \cr
#' }
#'
#' \code{mode_y = 'mesh'} therefore never reports \code{NA}: a face is always
#' testable. (The scan does skip faces flagged as deleted, but a surface is
#' rebuilt from its vertex and face matrices on every call, so no face arrives
#' flagged.)
#'
#' An empty chain cannot collide with anything, so reporting \code{FALSE} would
#' not be wrong in the set-theoretic sense. It is avoided for two reasons. In
#' \code{'points'} mode \code{hit_unit} is aligned to the rows of \code{y} so
#' that it can be indexed straight back into the input matrix, which means
#' separator rows occupy slots; calling those \code{FALSE} asserts something
#' about a point that does not exist. And \code{sum(is.na(hit_unit))} counts the
#' degenerate units that were passed in -- empty or single-vertex streamlines in
#' a tract file, deleted faces in a surface -- which is usually a data-quality
#' signal worth seeing rather than silently folding into "missed". Once
#' flattened to \code{FALSE} that count cannot be recovered.
#'
#' Most uses never have to deal with it. \code{which} already drops \code{NA},
#' so \code{which(hit_unit)} gives the colliding units directly, and
#' \code{representation$unit} is that same vector without any coercion.
#' \code{NA} only intrudes on \code{sum} / \code{any} / \code{if} without
#' \code{na.rm}, and on logical subsetting such as \code{tracts[hit_unit]},
#' which yields \code{NULL} entries. Where a plain mask really is wanted:
#'
#' \preformatted{
#'   mask <- res$hit_unit \%in\% TRUE     # NA becomes FALSE
#' }
#'
#' \code{test_level = 'whole'} reports \emph{a} colliding unit, not necessarily
#' the first one in order: the scan stops as soon as any thread finds a
#' collision. \code{collide} itself is deterministic, as are the two lower
#' levels, which measure each unit in row order.
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
#' # Which streamlines touch the sphere, and where does each come closest?
#' result$hit_unit
#' result$representation
#'
#' # What was actually asked, including what `mode = "auto"` resolved to
#' result$summary
#'
#' # A bundle assembled by appending a separator to every streamline gives one
#' # unit per streamline, even where a streamline is empty: here the middle one
#' # contributes only its own separator, and reports NA rather than vanishing
#' tracts <- list(
#'   cbind(seq(-3, 3, by = 0.5), 0, 0),
#'   matrix(numeric(0), ncol = 3),
#'   cbind(seq(-3, 3, by = 0.5), 5, 0)
#' )
#' bundle <- do.call("rbind", lapply(tracts, function(line) rbind(line, NA)))
#' vcg_detect_collision(roi, bundle, mode_y = "segments",
#'                      radius = 0.1)$hit_unit
#'
#' # The same question, giving up on a streamline once it is known to collide.
#' # `hit_unit` is unchanged; `index` names the first colliding segment instead
#' # of the closest one, which here happens to be the same segment
#' vcg_detect_collision(roi, streamlines, mode_y = "segments", radius = 0.1,
#'                      test_level = "unit")$representation
#'
#' # Does the bundle touch the sphere at all?
#' vcg_detect_collision(roi, streamlines, mode_y = "segments", radius = 0.1,
#'                      test_level = "whole")$collide
#'
#' # A point cloud against the same region, asking only for proximity
#' points <- rbind(c(0, 0, 0), c(1.05, 0, 0), c(10, 10, 10))
#' vcg_detect_collision(roi, points, radius = 0.1)$hit_unit
#'
#' # The centre of the sphere is far from its surface, so it only counts as a
#' # collision when the interior is included
#' vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1)$collide
#' vcg_detect_collision(roi, rbind(c(0, 0, 0)), radius = 0.1,
#'                      include_interior = TRUE)$collide
#'
#' @export
vcg_detect_collision <- function(
    x, y,
    mode_x = c("auto", "points", "segments", "mesh"),
    mode_y = c("auto", "points", "segments", "mesh"),
    radius = 0,
    test_level = c("element", "unit", "whole"),
    include_interior = FALSE) {

  mode_x <- match.arg(mode_x)
  mode_y <- match.arg(mode_y)
  test_level <- match.arg(test_level)

  radius <- as.double(radius)
  stopifnot2(
    length(radius) == 1L && !is.na(radius) && is.finite(radius) && radius >= 0,
    msg = "`vcg_detect_collision`: `radius` must be a single finite non-negative number."
  )

  include_interior <- isTRUE(as.logical(include_interior))

  side_x <- resolve_collision_input(x, mode_x, "x")
  side_y <- resolve_collision_input(y, mode_y, "y")

  if (include_interior && side_x$mode != "mesh") {
    stop("`vcg_detect_collision`: `include_interior = TRUE` needs `x` to be a closed surface, so `mode_x` must resolve to \"mesh\" (it resolved to \"", side_x$mode, "\").")
  }

  mode_code <- c(points = 0L, segments = 1L, mesh = 2L)
  level_code <- c(element = 0L, unit = 1L, whole = 2L)
  unit_type <- c(points = "point", segments = "line", mesh = "face")

  res <- vcgDetectCollision(
    x_ = side_x$vb,
    x_it_ = side_x$it,
    y_ = side_y$vb,
    y_it_ = side_y$it,
    mode_x = mode_code[[side_x$mode]],
    mode_y = mode_code[[side_y$mode]],
    radius = radius,
    test_level = level_code[[test_level]],
    include_interior = include_interior
  )

  # The C++ side hands back the hit table as four parallel vectors.
  representation <- data.frame(
    unit = res$unit,
    index = res$index,
    distance = res$distance,
    x_index = res$x_index
  )

  # Records the resolved modes, so a stored result can be read back without the
  # call that made it. The unit counts come from C++ because only it knows how
  # the separator rows actually grouped.
  summary <- list(
    x = list(
      mode = side_x$mode,
      unit_type = unit_type[[side_x$mode]],
      n_units = res$n_units_x
    ),
    y = list(
      mode = side_y$mode,
      unit_type = unit_type[[side_y$mode]],
      n_units = res$n_units_y
    ),
    radius = radius,
    test_level = test_level,
    include_interior = include_interior
  )

  if (test_level == "whole") {
    return(list(
      collide = res$collide,
      representation = representation,
      summary = summary
    ))
  }

  # `hit_unit` keeps its NAs on purpose; do not flatten them to FALSE here. NA
  # means "no geometry to test", FALSE means "tested and missed", and collapsing
  # the two would claim a separator row is a point that missed and would erase
  # the count of degenerate units. See the "Why hit_unit reports NA" section
  # above. A caller wanting a plain mask writes `res$hit_unit %in% TRUE`.
  list(
    collide = res$collide,
    hit_unit = res$hit_unit,
    representation = representation,
    summary = summary
  )
}
