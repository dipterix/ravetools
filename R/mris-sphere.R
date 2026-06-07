#' @title Project a surface onto a sphere and relax metric distortion
#' @description
#' Projects a closed triangular surface mesh radially onto a sphere and then
#' iteratively relaxes the metric distortion this introduces ("unfolding"), so
#' that geodesic distances and face orientations stay close to those of the
#' input surface - a step typically used to prepare an inflated cortical
#' surface for spherical registration.
#'
#' @details
#' \strong{This is a reduced procedure, not a complete reproduction of any
#' particular reference implementation.} The full unfolding procedure
#' described in the literature integrates seven weighted energy terms against
#' a separate reference surface through a multi-resolution, multi-thousand
#' line optimization pipeline, which is out of scope for this package. Instead,
#' this implementation keeps the two dominant terms and integrates them with
#' the same momentum-based machinery \code{\link{mris_inflate}} uses:
#' \describe{
#'   \item{Distance term (\code{l_dist})}{Restoring force pulling each vertex
#'     back towards the input mesh's own original distances to its
#'     neighbors.}
#'   \item{Area term (\code{l_area})}{Repulsive force that acts only on
#'     folded/negative-area faces, pushing their vertices apart - this is the
#'     actual "unfolding" mechanism.}
#' }
#' Both terms are integrated via momentum integration with gradient averaging
#' and a 1 mm per-step displacement cap, exactly as \code{\link{mris_inflate}}
#' does.
#'
#' A further simplification: the reference procedure relaxes a freshly
#' spherical-projected surface against a separate, previously-loaded
#' white-matter reference surface for the distance term; since this function
#' takes a single mesh as input, the input mesh's own metric (captured before
#' projection) is used as the reference instead - the same convention
#' \code{\link{mris_inflate}} already uses, and the practical analogue (the
#' white-matter surface is normally what gets inflated to produce a typical
#' input to this kind of unfolding step).
#'
#' @param mesh triangular mesh of class \code{'mesh3d'}, typically an inflated
#'   surface (see \code{\link{mris_inflate}}).  Must be closed, manifold, and
#'   genus-0 - the same hard precondition as \code{\link{mris_inflate}}.
#'   \code{mris_sphere} will raise an error if such defects are detected.
#' @param target_radius radius of the target sphere.  Default \code{100}.
#' @param n_averages number of gradient-averaging passes applied to the
#'   distance-term gradient each iteration (the same neighborhood-averaging
#'   mechanism \code{\link{mris_inflate}} uses).  Default \code{64}; reduced
#'   from the much larger averaging size used in production pipelines (tuned
#'   for cortical surfaces with hundreds of thousands of vertices) to a size
#'   more practical for typical \code{'mesh3d'} inputs.
#' @param niterations number of unfolding iterations.  Default \code{25}.
#' @param l_dist distance-preservation coefficient.  Default \code{1.0}.
#' @param l_area folded-face repulsion coefficient.  Default \code{1.0}.
#' @param momentum momentum coefficient.  Default \code{0.9}.
#' @param dt time step.  Default \code{0.05}.
#' @param verbose logical; print per-iteration progress.  Default \code{FALSE}.
#'
#' @returns The projected and relaxed surface as a \code{'mesh3d'} object with
#'   \code{vb}, \code{it}, and \code{normals}.
#'
#' @examples
#' if (is_not_cran()) {
#'
#' sphere <- vcg_sphere(sub_division = 3L)
#'
#' # deform so there is metric distortion to relax
#' sphere$vb[1, ] <- sphere$vb[1, ] * (1 + 0.2 * rnorm(ncol(sphere$vb)))
#'
#' # Fix defects
#' sphere <- vcg_fix_defects(sphere)
#'
#' result <- mris_sphere(
#'   sphere,
#'   n_averages = 8L,
#'   niterations = 10L,
#'   target_radius = 1,
#'   verbose = TRUE
#' )
#'
#' plot_mesh_polygon(list(sphere, result),
#'                   col = list("gray", "red"),
#'                   alpha = c(0.3, 0.9))
#'
#'
#'
#' }
#'
#' @export
mris_sphere <- function(
    mesh,
    target_radius = 100,
    n_averages    = 64L,
    niterations   = 25L,
    l_dist        = 1.0,
    l_area        = 1.0,
    momentum      = 0.9,
    dt            = 0.05,
    verbose       = FALSE
) {
    mesh <- meshintegrity(mesh, facecheck = TRUE)

    vb <- mesh$vb[1:3, , drop = FALSE]
    storage.mode(vb) <- "double"

    it <- mesh$it
    storage.mode(it) <- "integer"

    tmp <- mrisSphere(
        vb_           = vb,
        it_           = it,
        target_radius = as.double(target_radius)[[1L]],
        n_averages    = as.integer(n_averages)[[1L]],
        niterations   = as.integer(niterations)[[1L]],
        l_dist        = as.double(l_dist)[[1L]],
        l_area        = as.double(l_area)[[1L]],
        momentum      = as.double(momentum)[[1L]],
        dt            = as.double(dt)[[1L]],
        verbose       = as.logical(verbose)[[1L]]
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
