#' @title Localize white-matter and \verb{pial} surfaces from an intensity volume
#' @description
#' Starting from a single closed surface mesh (typically a smoothed estimate
#' of the white-matter boundary, in the same physical coordinate space as
#' \code{volume}) and a co-registered intensity volume (such as a normalized
#' \verb{T1} scan), iteratively deforms the surface in two passes to localize
#' the white/gray-matter and gray-matter/\verb{CSF} tissue-intensity
#' boundaries, producing a \code{white} and a \code{pial} surface.
#'
#' @details
#' The implementation keeps the two dominant ideas of the surface-placement
#' procedure described in the literature (see \strong{References}):
#' \enumerate{
#'   \item \strong{Intensity-target localization}: for each vertex, sample
#'     \code{volume} along the vertex's current normal at offsets spanning
#'     \eqn{\pm}\code{max_thickness} in steps of \code{step_size}, and pull
#'     the vertex toward the offset whose sampled intensity is closest to a
#'     single target value (\code{white_intensity} for the white-surface pass,
#'     \code{pial_intensity} for the \verb{pial}-surface pass).
#'   \item \strong{Smoothness}: a 1-ring Laplacian spring keeps the mesh
#'     regular while the per-vertex intensity term, which reacts independently
#'     to noisy image data, pulls vertices toward the tissue boundary.
#' }
#' Both terms are integrated with the same gradient-averaging and
#' momentum-integration machinery \code{\link{mris_inflate}} and
#' \code{\link{mris_sphere}} use: each inner iteration clears the gradient,
#' adds the intensity-target term, smooths it over \code{n_averages} passes of
#' 1-ring averaging, adds the locally-acting smoothness term (which is not
#' itself smoothed),
#' then takes a momentum-integration step with a 1 \verb{mm} per-step
#' displacement cap, and refreshes the vertex normals.
#'
#' The white-surface pass runs first, for \code{niterations} iterations; the
#' \verb{pial}-surface pass then continues from its result, with the momentum
#' velocity reset to rest, toward \code{pial_intensity} for another
#' \code{niterations} iterations.
#'
#' This is a \emph{reduced} port: the literature's procedure is a
#' multi-resolution optimization over roughly seven weighted energy terms
#' (intensity, intensity gradient, smoothness, self-intersection repulsion,
#' curvature, and more), using per-vertex gray/white/\verb{CSF} intensity
#' statistics derived from a prior segmentation. Reproducing that faithfully
#' is out of scope for this package; \code{white_intensity} and
#' \code{pial_intensity} are supplied directly here instead, for example the
#' midpoints between the typical white-matter/gray-matter and
#' gray-matter/\verb{CSF} intensities of \code{volume}.
#'
#' @param mesh triangular mesh of class \code{'mesh3d'} (or coercible via
#' \code{\link{ensure_mesh3d}}), in the same physical (\verb{RAS}) coordinate
#' space as \code{volume}. As with \code{\link{mris_inflate}}, the mesh
#' \strong{must} be closed, manifold, and genus-0 (watertight, no boundary or
#' non-manifold edges, single connected component); \code{mris_make_surfaces}
#' raises an error if such defects are detected
#' @param volume a 3-dimensional numeric array of image intensities (such as
#' a normalized \verb{T1} volume), co-registered with \code{mesh}
#' @param white_intensity target intensity for the white/gray-matter boundary;
#' the value the white-surface pass searches for along each vertex normal. For
#' a normalized \verb{T1} volume this is typically close to the midpoint
#' between the white-matter and gray-matter intensities
#' @param pial_intensity target intensity for the gray-matter/\verb{CSF}
#' boundary, analogous to \code{white_intensity} for the \verb{pial}-surface pass
#' @param IJK2RAS volume \verb{IJK} (zero-indexed voxel index) to surface
#' \verb{tkrRAS} transform (a \code{4x4} matrix); default is \code{NULL},
#' which assumes \code{volume} is a conformed volume (\verb{LIA} orientation,
#' \verb{1mm} isotropic 'voxels', centered at the volume midpoint, the same
#' convention \code{\link{fill_surface}} defaults to) and derives the
#' transform from \code{dim(volume)}; set this explicitly when \code{volume}
#' uses a different orientation or voxel size
#' @param max_thickness half-width, in \verb{mm}, of the per-vertex normal
#' search window the intensity-target term samples. Default \code{5}
#' @param step_size sampling step, in \verb{mm}, along the search window.
#' Default \code{0.4}
#' @param n_averages number of 1-ring gradient-averaging passes applied to the
#' intensity-target gradient before the smoothness term is added (mirrors
#' \code{\link{mris_inflate}}'s gradient-averaging step). Default \code{4}
#' @param niterations number of deformation iterations \emph{per pass}
#' (white, then \verb{pial}). Default \code{10}
#' @param l_intensity intensity-target term coefficient. Default \code{1.0}
#' @param l_spring smoothness (1-ring Laplacian spring) term coefficient.
#' Default \code{0.5}
#' @param momentum momentum coefficient. Default \code{0.9}
#' @param dt time step. Default \code{0.5}
#' @param verbose logical; print per-pass progress. Default \code{FALSE}
#'
#' @returns A named list of two \code{'mesh3d'} surfaces (each with
#' \code{vb}, \code{it}, and \code{normals}):
#' \describe{
#'   \item{\code{white}}{Surface localized to \code{white_intensity}.}
#'   \item{\code{pial}}{Surface localized to \code{pial_intensity}, continuing
#'     from \code{white}.}
#' }
#'
#' @references
#' Cortical surface-based analysis I: Segmentation and surface
#' reconstruction. \emph{NeuroImage}, 9(2), 179-194 (1999).
#'
#' @examples
#'
#' if (is_not_cran()) {
#'
#'
#' data("left_hippocampus_mask")
#' n_vox <- length(left_hippocampus_mask)
#'
#' volume <- left_hippocampus_mask + runif(n = n_vox, 0, 1)
#' vox2ras <- diag(1, 4)
#' mesh <- vcg_isosurface(volume, threshold_lb = 0.99)
#'
#' plot(mesh)
#'
#' # Fix defects
#' mesh <- vcg_fix_defects(mesh, verbose = TRUE, merge_tolerance = 1.75)
#'
#'
#' res <- mris_make_surfaces(
#'   mesh,
#'   volume,
#'   pial_intensity = 1.1,
#'   white_intensity = 1,
#'   IJK2RAS = vox2ras
#' )
#'
#' plot(res$pial)
#'
#'
#' }
#'
#'
#' @export
mris_make_surfaces <- function(
    mesh,
    volume,
    white_intensity,
    pial_intensity,
    IJK2RAS       = NULL,
    max_thickness = 5,
    step_size     = 0.4,
    n_averages    = 4L,
    niterations   = 10L,
    l_intensity   = 1.0,
    l_spring      = 0.5,
    momentum      = 0.9,
    dt            = 0.5,
    verbose       = FALSE
) {
  mesh <- meshintegrity(mesh, facecheck = TRUE)

  vb <- mesh$vb[1:3, , drop = FALSE]
  storage.mode(vb) <- "double"

  it <- mesh$it
  storage.mode(it) <- "integer"

  volume <- as.array(volume)
  if (length(dim(volume)) != 3L) {
    stop("`mris_make_surfaces`: `volume` must be a 3-dimensional array")
  }
  storage.mode(volume) <- "double"

  if (length(IJK2RAS) != 16L) {
    # Conformed-volume ('LIA', 1mm isotropic, centered) default; the same
    # convention `fill_surface` assumes when `IJK2RAS` is not supplied,
    # generalized here from a fixed cube edge to `dim(volume)`.
    vol_dim <- dim(volume)
    IJK2RAS <- matrix(
      nrow = 4, byrow = TRUE,
      c(-1, 0, 0, vol_dim[[1]] / 2,
        0, 0, 1, -vol_dim[[3]] / 2,
        0, -1, 0, vol_dim[[2]] / 2,
        0, 0, 0, 1))
  } else {
    IJK2RAS <- matrix(as.double(IJK2RAS), nrow = 4)
  }

  ras2ijk <- solve(IJK2RAS)
  storage.mode(ras2ijk) <- "double"

  tmp <- mrisMakeSurfaces(
    vb_             = vb,
    it_             = it,
    volume_         = volume,
    ras2ijk_        = ras2ijk,
    white_intensity = as.double(white_intensity)[[1L]],
    pial_intensity  = as.double(pial_intensity)[[1L]],
    max_thickness   = as.double(max_thickness)[[1L]],
    step_size       = as.double(step_size)[[1L]],
    n_averages      = as.integer(n_averages)[[1L]],
    niterations     = as.integer(niterations)[[1L]],
    l_intensity     = as.double(l_intensity)[[1L]],
    l_spring        = as.double(l_spring)[[1L]],
    momentum        = as.double(momentum)[[1L]],
    dt              = as.double(dt)[[1L]],
    verbose         = as.logical(verbose)[[1L]]
  )

  pack <- function(part) {
    structure(
      list(
        vb      = rbind(part$vb, 1),
        it      = it,
        normals = rbind(part$normals, 1)
      ),
      class = c("ravetools_mesh3d", "mesh3d")
    )
  }

  list(
    white = pack(tmp$white),
    pial  = pack(tmp$pial)
  )
}
