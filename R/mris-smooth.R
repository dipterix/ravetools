#' @title Smooth a surface mesh
#' @description
#' Smooths a closed triangular surface mesh by repeatedly replacing each
#' vertex with the average of itself and its immediate neighbors
#' (Laplacian/neighbor-averaging smoothing), optionally rescaling the surface
#' back to its original area afterwards.
#'
#' @details
#' For each of \code{npasses} passes, the algorithm:
#' \enumerate{
#'   \item Replaces each vertex position by the include-self mean of itself
#'     and its 1-ring (directly-connected) neighbors, repeated
#'     \code{niterations} times.
#'   \item Recomputes vertex normals and total surface area.
#'   \item Optionally rescales the surface back to its original area (only if
#'     \code{rescale = TRUE}).
#' }
#'
#' @param mesh triangular mesh of class \code{'mesh3d'}.  Must be closed,
#'   manifold, and genus-0 (watertight, no boundary or non-manifold edges,
#'   single connected component) - the same hard precondition as
#'   \code{\link{mris_inflate}} (and for the same reason: an open or
#'   defective mesh introduces a persistent directional bias into the
#'   1-ring averaging).  \code{mris_smooth} will raise an error if such
#'   defects are detected.
#' @param niterations number of include-self 1-ring averaging rounds applied
#'   to the vertex positions in each pass.  Default \code{10}.
#' @param npasses number of outer smoothing passes.  Default \code{1}.
#' @param rescale logical; whether to rescale the surface back to its original
#'   area after each pass.  Default \code{FALSE}.
#' @param verbose logical; print per-pass progress.  Default \code{FALSE}.
#'
#' @returns The smoothed surface as a \code{'mesh3d'} object with \code{vb},
#'   \code{it}, and \code{normals}.
#'
#' @examples
#' if (is_not_cran()) {
#'
#' sphere <- vcg_sphere(sub_division = 4L)
#'
#' # roughen the sphere slightly so smoothing has something to do
#' sphere$vb[1, ] <- sphere$vb[1, ] * (1 + 0.05 * rnorm(ncol(sphere$vb)))
#'
#' smoothed <- mris_smooth(sphere, niterations = 5L, verbose = TRUE)
#'
#' plot_mesh_polygon(
#'   list(sphere, smoothed),
#'   alpha = c(0.3, 0.5),
#'   col = list("gray", "red"),
#'   main = "Gray: original; red: smoothed"
#' )
#'
#' }
#'
#' @export
mris_smooth <- function(
    mesh,
    niterations = 10L,
    npasses     = 1L,
    rescale     = FALSE,
    verbose     = FALSE
) {
    mesh <- meshintegrity(mesh, facecheck = TRUE)

    vb <- mesh$vb[1:3, , drop = FALSE]
    storage.mode(vb) <- "double"

    it <- mesh$it
    storage.mode(it) <- "integer"

    tmp <- mrisSmooth(
        vb_         = vb,
        it_         = it,
        niterations = as.integer(niterations)[[1L]],
        npasses     = as.integer(npasses)[[1L]],
        rescale     = as.logical(rescale)[[1L]],
        verbose     = as.logical(verbose)[[1L]]
    )

    structure(
        list(
            vb      = rbind(tmp$vb, 1),
            it      = it,
            normals = rbind(tmp$normals, 1)
        ),
        class = "mesh3d"
    )
}
