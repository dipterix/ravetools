#' @title Isotropic re-triangulation of a surface mesh
#' @description
#' Re-triangulates a closed surface mesh so that all edges approach a uniform
#' target length, improving triangle regularity and quality without changing
#' the surface topology or shape. The output has a different vertex and face
#' count than the input but closely tracks the original surface geometry.
#'
#' @details
#' Each iteration applies three steps following the procedure described in the
#' \strong{References}:
#' \enumerate{
#'   \item \strong{Edge split}: every edge longer than
#'     \eqn{4/3 \times \code{target\_edge\_length}} is split at its midpoint.
#'     One-edge splits produce 2 sub-triangles; two-edge splits produce 3;
#'     three-edge splits produce 4 (the standard 1-to-4 uniform refinement).
#'   \item \strong{Edge collapse}: every edge shorter than
#'     \eqn{4/5 \times \code{target\_edge\_length}} is collapsed to its
#'     midpoint. A collapse is skipped if it would flip any surrounding
#'     face normal (manifold-safety check).
#'   \item \strong{Tangential smoothing}: \code{n_smooth} passes of
#'     uniform-weight 1-ring \verb{Laplacian} averaging, with each
#'     displacement projected onto the vertex tangent plane (the normal
#'     component is removed) before application. This improves vertex
#'     regularity while keeping vertices near the original surface.
#' }
#' After all iterations, vertex normals are refreshed.
#'
#' Unlike \code{\link{vcg_uniform_remesh}}, which uses volumetric resampling
#' and may change topology, this function operates entirely on the mesh
#' surface and preserves the input genus and manifold structure.
#'
#' @param mesh triangular mesh of class \code{'mesh3d'} (or coercible via
#' \code{\link{ensure_mesh3d}}). The mesh \strong{must} be closed, manifold,
#' and genus-0 (watertight, no boundary or non-manifold edges, single
#' connected component); \code{mris_remesh} raises an error if such defects
#' are detected
#' @param target_edge_length numeric; desired uniform edge length in the same
#' units as \code{mesh$vb}. Default \code{NULL} uses the current average edge
#' length (quality improvement without size change)
#' @param niterations integer; number of split/collapse/smooth iterations.
#' Default \code{5}
#' @param n_smooth integer; number of tangential-smooth passes per iteration.
#' Default \code{2}
#' @param damping numeric; fraction of the tangential displacement applied at
#' each smooth pass (between 0 and 1). Default \code{0.99}
#' @param verbose logical; print per-iteration vertex and face counts.
#' Default \code{FALSE}
#'
#' @returns A \code{'mesh3d'} object with \code{vb}, \code{it}, and
#' \code{normals}. Vertex and face counts differ from the input; the surface
#' geometry is preserved.
#'
#' @references
#' A \verb{remeshing} approach to \verb{multiresolution} modeling.
#' \emph{Proceedings of Shape \verb{Modelling} International}, 49-58 (2003).
#'
#'
#'
#' sphere <- vcg_sphere()
#' sphere
#'
#' vcg_average_edge_length(sphere)
#' plot(sphere)
#'
#' remeshed <- mris_remesh(
#'   sphere,
#'   target_edge_length = 0.3
#' )
#'
#' plot(remeshed)
#'
#'
#'
#' @export
mris_remesh <- function(
    mesh,
    target_edge_length = NULL,
    niterations        = 5L,
    n_smooth           = 2L,
    damping            = 0.99,
    verbose            = FALSE
) {
  mesh <- meshintegrity(mesh, facecheck = TRUE)

  vb <- mesh$vb[1:3, , drop = FALSE]
  storage.mode(vb) <- "double"

  it <- mesh$it
  storage.mode(it) <- "integer"

  if (is.null(target_edge_length)) {
    target_edge_length <- vcg_average_edge_length(mesh)
  }

  tmp <- mrisRemesh(
    vb_            = vb,
    it_            = it,
    target_edge_len = as.double(target_edge_length)[[1L]],
    niterations    = as.integer(niterations)[[1L]],
    n_smooth       = as.integer(n_smooth)[[1L]],
    damping        = as.double(damping)[[1L]],
    verbose        = as.logical(verbose)[[1L]]
  )

  structure(
    list(
      vb      = rbind(tmp$vb, 1),
      it      = tmp$it,
      normals = rbind(tmp$normals, 1)
    ),
    class = c("ravetools_mesh3d", "mesh3d")
  )
}
