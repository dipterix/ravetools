#' @title Convex hull of a point cloud or surface
#' @description
#' Computes the smallest convex closed surface that contains a set of points
#' (for a surface, its vertices), with the \verb{Quickhull} implementation of
#' the bundled \verb{VCG} library or with a port of the \verb{quickhull3d}
#' library to 'C++'.
#' @param x a matrix of point coordinates with 3 columns and one row per
#' point, or a surface (see 'Coercing Surface Inputs' below), in which case
#' its vertices are used and its faces ignored
#' @param method algorithm: \code{"vcg"} (default) for the \verb{Quickhull}
#' of \verb{VCG} (\verb{Barber} and colleagues, 1996), or
#' \code{"quickhull3d"} for a port of the \verb{quickhull3d} library, which
#' follows \verb{Lloyd}'s \verb{QuickHull3D} and merges neighboring faces
#' that lie in nearly the same plane; see 'Details'
#' @param tolerance non-negative distance, in the units of \code{x}: a point
#' closer than this to the plane of a hull face counts as lying on the face
#' rather than outside it. The default \code{NA} uses the method's own value:
#' \code{1e-10} times the diagonal of the bounding box for \code{"vcg"}, and
#' for \code{"quickhull3d"} the rounding error of the coordinates, three
#' times the machine precision times the sum, over the three axes, of the
#' largest absolute coordinate. Values below that rounding error are raised
#' to it, because below it the visibility tests contradict each other
#' @returns A closed triangular mesh of class \code{'mesh3d'} with:
#' \item{\code{vb}}{coordinates of the hull vertices, each one exactly an
#' input point, in the input order}
#' \item{\code{it}}{triangular faces, ordered counter-clockwise when seen
#' from outside (face normals point outward)}
#' \item{\code{normals}}{area-weighted vertex normal vectors}
#' and attribute \code{"vertex_index"}: for each hull vertex, the row of
#' \code{x} (or the vertex index of the surface) it comes from, increasing;
#' \code{attr(hull, "vertex_index")[hull$it]} indexes the faces by input
#' point.
#'
#' @details
#' Both methods start from a tetrahedron of extreme points and repeatedly add
#' the point farthest outside a face, replacing the faces that point can see.
#' Points within \code{tolerance} of a face plane count as on the face.
#'
#' The \verb{VCG} implementation keeps every face a triangle, so a flat region
#' of the hull (common for points on a voxel grid) is split into arbitrary
#' triangles, and it relies on the tolerance alone to stay consistent under
#' rounding. The \verb{quickhull3d} method also merges neighboring faces that
#' lie in one plane or are not convex within the tolerance into polygons,
#' which it splits into triangles fanning out from one corner; such a fan can
#' contain triangles of zero area when a polygon has three corners on a line.
#'
#' Points that are identical, or all lie on one line or in one plane, enclose
#' no volume and raise an error. So does a result that fails the check for a
#' closed, consistently oriented surface, a numerical failure that a larger
#' \code{tolerance} or the other method can avoid.
#'
#' @references
#' \verb{Barber} CB, \verb{Dobkin} DP, \verb{Huhdanpaa} H (1996). The
#' \verb{quickhull} algorithm for convex hulls. \emph{\verb{ACM} Transactions on
#' Mathematical Software}, 22(4), 469-483.
#'
#' \verb{Lloyd JE}. \verb{QuickHull3D}: a robust 3D convex hull algorithm in
#' Java. \url{https://www.cs.ubc.ca/~lloyd/java/quickhull3d.html}
#'
#' \verb{Poppe} M. \verb{quickhull3d}: a robust \verb{quickhull}
#' implementation to find the convex hull of a set of 3D points.
#' \url{https://github.com/mauriciopoppe/quickhull3d}
#'
#' @inheritSection ensure_mesh3d Coercing Surface Inputs
#'
#' @examples
#'
#' # the 8 corners of a unit cube and 100 points inside it
#' set.seed(1)
#' x <- rbind(
#'   as.matrix(expand.grid(0:1, 0:1, 0:1)),
#'   matrix(runif(300), ncol = 3)
#' )
#' hull <- vcg_convex_hull(x)
#'
#' # the rows of `x` on the hull are the corners
#' attr(hull, "vertex_index")
#'
#' vcg_mesh_volume(hull)
#'
#' if(is_not_cran()) {
#'
#' data("left_hippocampus_mask")
#' mesh <- vcg_isosurface(left_hippocampus_mask)
#' hull <- vcg_convex_hull(mesh, method = "quickhull3d")
#'
#' rgl_view({
#'   rgl_call("shade3d", mesh, col = 2)
#'   rgl_call("wire3d", hull, col = 3)
#' })
#'
#' }
#'
#' @export
vcg_convex_hull <- function(x, method = c("vcg", "quickhull3d"),
                            tolerance = NA) {
  method <- match.arg(method)
  xyz <- as_point_cloud_matrix(x, "vcg_convex_hull")
  if (ncol(xyz) < 4L) {
    stop("vcg_convex_hull: at least 4 points are needed for a 3D convex hull")
  }
  if (!all(is.finite(xyz))) {
    stop("vcg_convex_hull: point coordinates must be finite (no `NA`, `NaN`, or `Inf`)")
  }
  # Both backends start from a tetrahedron of extreme points, which stays flat
  # when the points enclose no volume (vcglib then reads uninitialized
  # pointers). The singular values of the centered points measure their
  # spread along the three principal axes.
  spread <- svd(xyz - rowMeans(xyz), nu = 0L, nv = 0L)$d
  if (spread[[3]] <= sqrt(.Machine$double.eps) * spread[[1]]) {
    stop("vcg_convex_hull: the points are identical, or all lie on one line or in one plane, so they enclose no volume")
  }
  tolerance <- as.double(tolerance)
  if (length(tolerance) != 1L || isTRUE(tolerance < 0)) {
    stop("vcg_convex_hull: `tolerance` must be `NA` or a non-negative number")
  }

  # Below the rounding error of the coordinates, the visibility tests of both
  # methods contradict each other on points that share planes (vcglib then
  # corrupted its face links and crashed; quickhull3d's merges never ended),
  # so smaller tolerances are raised to it. quickhull3dConvexHull does that
  # itself, with quickhull3d's own computation of the same bound.
  rounding <- 3 * .Machine$double.eps * sum(apply(abs(xyz), 1L, max))

  faces <- tryCatch(
    switch(
      method,
      "vcg" = {
        if (is.na(tolerance)) {
          extent <- apply(xyz, 1L, function(v) { diff(range(v)) })
          tolerance <- 1e-10 * sqrt(sum(extent^2))
        }
        vcgConvexHull(xyz, max(tolerance, rounding))
      },
      "quickhull3d" = {
        # a negative value makes the port use quickhull3d's own tolerance
        if (is.na(tolerance)) { tolerance <- -1 }
        quickhull3dConvexHull(xyz, tolerance)
      }
    ),
    error = function(e) {
      stop(sprintf(paste(
        "vcg_convex_hull: method '%s' failed (%s); try a larger `tolerance`",
        "or the other `method`"), method, conditionMessage(e)), call. = FALSE)
    }
  )

  vertex_index <- sort(unique(as.vector(faces)))
  it <- matrix(match(faces, vertex_index), nrow = 3L)
  if (!is_closed_surface(it, length(vertex_index))) {
    stop(sprintf(paste(
      "vcg_convex_hull: method '%s' produced a surface that is not closed",
      "and consistently oriented (a numerical failure); try a larger",
      "`tolerance` or the other `method`"), method))
  }

  hull <- structure(
    list(
      vb = rbind(xyz[, vertex_index, drop = FALSE], 1, deparse.level = 0),
      it = it
    ),
    class = c("ravetools_mesh3d", "mesh3d")
  )
  # area weighting: the zero-area triangles a merged polygon's fan can contain
  # add nothing, whereas angle weighting would give them a weight up to pi and
  # an arbitrary normal
  hull <- vcg_update_normals(hull, weight = "area")
  attr(hull, "vertex_index") <- vertex_index
  hull
}

# Whether triangles `it` (vertex indices 1 to `n_vertex`, all used) form a
# closed, consistently oriented surface of genus 0: every edge borders two
# faces that traverse it in opposite directions, and V - E + F = 2
is_closed_surface <- function(it, n_vertex) {
  from <- c(it[1L, ], it[2L, ], it[3L, ])
  to <- c(it[2L, ], it[3L, ], it[1L, ])
  directed <- from * (n_vertex + 1) + to
  if (anyDuplicated(directed)) { return(FALSE) }
  if (!all((to * (n_vertex + 1) + from) %in% directed)) { return(FALSE) }
  n_vertex - length(directed) / 2 + ncol(it) == 2
}
