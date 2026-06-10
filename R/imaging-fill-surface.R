#' @title Fill a volume cube based on water-tight surface
#' @author Zhengjia Wang
#' @description Create a cube volume (\code{256} 'voxels' on each margin), fill
#' in the 'voxels' that are inside of the surface.
#'
#' @param surface a surface mesh; accepted classes include \code{'mesh3d'}
#' (from \code{rgl}), \code{'fs.surface'} (from \code{freesurferformats}),
#' \code{'ieegio_surface'} (from \pkg{ieegio}), or a bare list containing a
#' \code{vb} vertex matrix; see \code{\link{ensure_mesh3d}} for details and
#' for how \code{'ieegio_surface'} inputs are coerced
#' @param inflate amount of 'voxels' to inflate on the final result; must be
#' a non-negative integer. A zero \code{inflate} value means the resulting
#' volume is tightly close to the surface
#' @param close_radius radius (in 'voxels' along the finest 'IJK' axis; see
#' \code{IJK2RAS}) of the morphological closing (dilate-then-erode) operation
#' used to connect the embedded surface into a water-tight shell \emph{and} to
#' remove small topological handles/tunnels (for example a hole drilled
#' through an otherwise solid structure) before the surface is filled; default
#' is \code{NULL}, which derives a radius from the surface geometry (capped at
#' \code{15} 'voxels'). Handles/tunnels wider than roughly twice this radius
#' (in physical, \verb{RAS} space) will not be closed, so when the input
#' surface is known to contain larger topological defects, set this to a
#' larger number explicitly. Larger radii cost more time and memory (the
#' closing operates on a dense \code{resolution^3} volume) and - because the
#' operation is performed in a fixed-size volume - must leave enough margin
#' that the dilated shell does not reach the volume boundary; radii that are
#' too large for the working cube will be capped with a warning
#' @param resolution number of 'voxels' along each margin of the working
#' cube volume that the surface is embedded into; default is \code{256}.
#' Larger values increase the spatial precision of the fill and the closing
#' operation (at the cost of memory and computation time, which scale with
#' \code{resolution^3}); the value must leave enough room for the surface
#' (transformed into 'IJK' space via \code{IJK2RAS}) plus the closing margin
#' to fit inside the cube
#' @param IJK2RAS volume 'IJK' (zero-indexed coordinate index) to
#' \code{'tkrRAS'} transform, default is automatically determined
#' leave it `NULL` if you don't know how to set it. When the linear part of
#' \code{IJK2RAS} encodes voxel sizes that differ across axes (\code{pixdim},
#' i.e. the three 'IJK' axes do not advance by the same physical distance per
#' 'voxel'), the dilate/erode operations underlying \code{close_radius} and
#' \code{inflate} automatically use a different number of 'voxels' along each
#' axis so that the resulting structuring element stays isotropic in physical
#' (\verb{RAS}) space
#' @param preview whether to preview the results; default is false
#' @param preview_frame integer from 1 to 256 the depth frame used to generate
#' preview.
#'
#' @details This function creates a volume (256 on each margin) and fill in
#' the volume from a surface mesh. The surface vertex points will be embedded
#' into the volume first. These points may not be connected together, hence
#' for each 'voxel', a cube patch will be applied to grow the volume. Then,
#' the volume will be bucket-filled from a corner, forming a negated mask of
#' "outside-of-surface" area. The inverted bucket-filled volume is then shrunk
#' so the mask boundary tightly fits the surface
#'
#'
#' @returns A list containing the filled volume and parameters used to generate
#' the volume
#'
#' @inheritSection ensure_mesh3d Coercing Surface Inputs
#'
#' @seealso \code{\link{ensure_mesh3d}}
#'
#' @examples
#'
#' \donttest{
#'
#' # takes > 5s to run example
#'
#' # Generate a sphere
#' surface <- vcg_sphere()
#' surface$vb[1:3, ] <- surface$vb[1:3, ] * 50
#'
#' fill_surface(surface, preview = TRUE)
#'
#' }
#'
#' @name fill_surface
NULL

#' @title Coerce a surface object to a \code{'mesh3d'} mesh
#' @description
#' Internal helper used throughout \pkg{ravetools} to accept a variety of
#' surface representations and return a list of class \code{'mesh3d'} (the
#' format used by the \code{rgl} package). Most mesh-consuming functions in
#' this package call \code{ensure_mesh3d} on their surface arguments, so the
#' coercion rules described here apply to all of them.
#'
#' @param surface a surface object. One of the following:
#' \describe{
#'   \item{\code{'mesh3d'}}{returned unchanged.}
#'   \item{\code{'fs.surface'}}{a surface in the format produced by
#'     \code{freesurferformats::read.fs.surface}. Vertices and faces are
#'     copied into a \code{'mesh3d'} list; zero-indexed faces are bumped
#'     by 1.}
#'   \item{\code{'ieegio_surface'}}{a surface object produced by
#'     \pkg{ieegio} (e.g. from \code{read_surface} or
#'     \code{volume_to_surface}). See the section below for important
#'     behavior changes.}
#'   \item{a bare \code{list}}{must contain a numeric \code{vb} matrix with
#'     3 or 4 rows; the list is reclassified as \code{'mesh3d'} with no
#'     other modification.}
#' }
#' Any other input triggers an error.
#'
#' @returns An object of class \code{'mesh3d'} with at least the \code{vb}
#' (vertex) component and, when face information is available, an \code{it}
#' (triangle index) component.
#'
#' @section Coercing Surface Inputs:
#' The surface objects are converted to \code{'mesh3d'} object before
#' applying further calculations.
#'
#' When \code{surface} is a surface \pkg{ieegio} object, the returned
#' \code{mesh3d$vb} contains vertices that have been left-multiplied by
#' \code{surface$geometry$transforms[[1]]} (the first transform stored in the
#' geometry, typically the \code{ScannerAnat} or voxel-to-world transform).
#'
#' \strong{Breaking change:} Earlier versions (before 0.2.6) of \pkg{ravetools}
#' returned the raw \code{surface$geometry$vertices} without applying any
#' transform, so downstream code often multiplied by
#' \code{surface$geometry$transforms[[1]]} (or an equivalent) manually before
#' working in world space. Such code will now \emph{double}
#' apply the transform and produce incorrect coordinates. If you previously
#' applied a transform from \code{surface$geometry$transforms} by hand after
#' calling a \pkg{ravetools} mesh function on an \code{'ieegio_surface'},
#' remove that manual step.
#'
#' Surfaces with an empty or missing \code{geometry$transforms} list (for
#' example, surfaces produced by \pkg{ieegio}'s \code{volume_to_surface},
#' which stores an identity transform) are unaffected.
#'
#' If \code{geometry$transforms} contains multiple transforms targeting
#' different coordinate spaces, only the first one is used. Callers that need
#' a specific target space should select and apply that transform themselves
#' before calling \pkg{ravetools} mesh functions.
#'
#' @examples
#'
#' # mesh3d input is returned unchanged in shape
#' sphere <- vcg_sphere()
#' m <- ensure_mesh3d(sphere)
#' identical(m$vb, sphere$vb)
#'
#' # A bare list with a `vb` slot is reclassified to mesh3d
#' bare <- list(vb = rbind(matrix(rnorm(30), nrow = 3), 1))
#' m <- ensure_mesh3d(bare)
#' inherits(m, "mesh3d")
#'
#' @export
ensure_mesh3d <- function(surface) {
  # check if surface is read from `freesurferformats`
  if (inherits(surface, "fs.surface")) {
    # surface <- freesurferformats::fs.surface.to.tmesh3d( surface )
    face <- t(surface$faces)
    if ( min(face, na.rm = TRUE) == 0 ) {
      face <- face + 1
    }
    surface <- structure(list(
      vb = rbind(t(surface$vertices), 1),
      it = face), class = c("ravetools_mesh3d", "mesh3d"))
  } else if (inherits(surface, "ieegio_surface")) {
    if (!inherits(surface, "ieegio_surface_contains_geometry")) {
      stop("Surface object has no geometry")
    }
    if (isTRUE(surface$sparse) || !length(surface$sparse_node_index)) {
      vertices <- surface$geometry$vertices
    } else {
      node_index <- surface$sparse_node_index
      vertices <- array(0.0, c(3, max(node_index)))
      vertices[1:3, ] <- surface$geometry$vertices[1:3, ]
      vertices <- rbind(vertices, 1)
    }

    # Apply transforms
    if (length(surface$geometry$transforms)) {
      transform <- surface$geometry$transforms[[1]]
      vertices <- transform %*% vertices
    }

    face <- surface$geometry$faces
    if (length(face)) {
      face_start <- surface$geometry$face_start
      if (length(face_start) == 1 && !isTRUE(face_start == 1)) {
        face <- face - face_start + 1L
      }
      storage.mode(face) <- "integer"

      surface <- structure(list(
        vb = vertices,
        it = face), class = c("ravetools_mesh3d", "mesh3d"))
    } else {
      surface <- structure(list(
        vb = vertices), class = c("ravetools_mesh3d", "mesh3d"))
    }
  }
  if (!inherits(surface, "mesh3d")) {
    # check if surface has names
    if (
      is.list(surface) &&
      "vb" %in% names(surface) &&
      is.matrix(surface$vb) &&
      nrow(surface$vb) %in% c(3L, 4L)
    ) {
      surface <- structure(surface, class = c("ravetools_mesh3d", "mesh3d"))
    } else {
      stop("ravetools:::ensure_mesh3d: `surface` must be a mesh3d object (package rgl) or a fs.surface object (package freesurferformats)")
    }
  }
  if (!inherits(surface, "ravetools_mesh3d")) {
    class(surface) <- c("ravetools_mesh3d", class(surface))
  }
  return(surface)
}

#' @rdname fill_surface
#' @export
fill_surface <- function(surface, inflate = 0, close_radius = NULL,
                         resolution = 256L, IJK2RAS = NULL,
                         preview = FALSE, preview_frame = 128) {
  # DIPSAUS DEBUG START
  # surface <- freesurferformats::read.fs.surface('~/Dropbox (PENN Neurotrauma)/RAVE/Samples/raw/PAV010/rave-imaging/fs/surf/rh.pial')
  # IJK2RAS <- NULL
  # verbose <- TRUE
  # preview <- TRUE
  # output_format <- "rgl"
  # inflate <- 3

  inflate <- as.integer(inflate)
  if (inflate < 0) { inflate <- 0 }

  # setup
  preview2D <- function(expr) {
    expr <- substitute(expr)
    if (interactive() && preview) {
      eval(expr)
    }
  }

  # Initialize the parameters
  resolution <- as.integer(resolution)
  if (length(resolution) != 1L || is.na(resolution) || resolution < 16L) {
    stop("`fill_surface`: `resolution` must be a single integer no less than 16")
  }
  if (length(IJK2RAS) != 16L) {
    IJK2RAS <- matrix(
      nrow = 4, byrow = TRUE,
      c(-1, 0, 0, resolution / 2,
        0, 0, 1, -resolution / 2,
        0, -1, 0, resolution / 2,
        0, 0, 0, 1))
  }

  # Voxel size ('pixdim') along each 'IJK' axis (in the units of `IJK2RAS`,
  # typically millimeters): the magnitude of each column of the linear part of
  # `IJK2RAS` is the physical distance spanned by one 'voxel' step along that
  # axis. The default `IJK2RAS` above is an axis-permuted identity (isotropic,
  # pixdim = c(1, 1, 1)); user-supplied transforms (e.g. from non-conformed
  # volumes) may be anisotropic.
  pixdim <- sqrt(colSums(IJK2RAS[1:3, 1:3, drop = FALSE] ^ 2))

  # Convert an isotropic radius in physical ('RAS') space - expressed as a
  # 'voxel' count along the finest (smallest-pixdim) axis - into a per-axis
  # 'voxel' count, so `grow_volume`'s box structuring element ends up
  # isotropic in physical space even when `pixdim` is anisotropic: an axis
  # with larger 'voxels' needs proportionally fewer of them to span the same
  # physical distance.
  to_voxel_radius <- function(radius) {
    as.integer(round(radius * min(pixdim) / pixdim))
  }

  # Make sure the surface is a `mesh3d` object
  surface_orig <- ensure_mesh3d(surface)

  # Transform surface to IJK voxel space so can be fitted into a volume
  # Remember IJK starts from 0
  surface <- surface_orig
  surface$vb <- solve(IJK2RAS) %*% surface$vb

  # Largest per-axis radius the dilated shell can use without reaching the
  # volume boundary (the bucket-fill corner must remain outside the dilated
  # shell)
  vertex_range <- apply(surface$vb, 1, range)
  safe_radius_per_axis <- floor(pmin(
    vertex_range[1, 1:3],
    resolution - vertex_range[2, 1:3]
  ) - 2)
  safe_radius <- min(safe_radius_per_axis)

  if (length(close_radius)) {
    close_radius <- as.numeric(close_radius)
    if (length(close_radius) != 1L || is.na(close_radius) || close_radius < 0) {
      stop("`fill_surface`: `close_radius` must be a single non-negative number")
    }
  } else {
    # Previous default behavior: derive a radius from the surface geometry,
    # capped at 15 'voxels' (along the finest axis)
    close_radius <- min(safe_radius, 15)
  }

  close_radius_voxels <- to_voxel_radius(close_radius)
  if (any(close_radius_voxels > safe_radius_per_axis)) {
    close_radius_voxels <- pmin(close_radius_voxels, safe_radius_per_axis)
    warning(sprintf(
      paste("`fill_surface`: requested `close_radius` (%g) is too large for",
            "this surface's position in the working %d^3 volume along at",
            "least one axis and has been capped at c(%s) 'voxels' (per I/J/K",
            "axis) to avoid the dilated shell reaching the volume boundary;",
            "topological handles/tunnels wider than roughly twice this",
            "radius (in physical space) may not be closed"),
      close_radius, resolution, paste(close_radius_voxels, collapse = ", ")
    ))
  }

  # Creating a volume
  volume <- array(0L, dim = rep(resolution, 3))

  # embed the surface in volume space
  surface_index <- round(surface$vb[c(1, 2, 3), ])
  surface_index <- surface_index[1, ] + surface_index[2, ] * resolution +
    surface_index[3, ] * (resolution^2) + 1
  volume[surface_index] <- 1L

  # Morphological closing (dilate by `close_radius`, then erode back): connects
  # the embedded surface points into a water-tight shell *and* bridges/closes
  # any topological handles or tunnels no wider than ~2x `close_radius` (in
  # physical space). The structuring element is `close_radius_voxels`, a
  # per-axis 'voxel' count that compensates for anisotropic `pixdim`.
  volume_grew <- grow_volume(
    volume,
    x = close_radius_voxels[[1]],
    y = close_radius_voxels[[2]],
    z = close_radius_voxels[[3]]
  )

  # bucket-fill from the corner, bucketFillVolume is in-place
  volume_filled <- bucketFillVolume(volume_grew + 0L, x = 1, y = 1, z = 1, fill = 1L)

  # Fill the voxels within the surface and shrink back by `close_radius`,
  # leaving a net `inflate` (also pixdim-corrected) relative to the surface
  shrink_amount_voxels <- close_radius_voxels - to_voxel_radius(inflate)
  volume2 <- 1L - grow_volume(
    volume_filled - volume_grew,
    x = shrink_amount_voxels[[1]],
    y = shrink_amount_voxels[[2]],
    z = shrink_amount_voxels[[3]]
  )

  # preview
  preview2D({
    oldPar <- graphics::par(c("mfrow", "mar"))
    graphics::par(mfrow = c(1, 3),
                  mar = c(3.1, 0.1, 3.1, 0.1))
    on.exit({ do.call(graphics::par, oldPar) })

    frame <- as.integer(preview_frame)
    if (is.na(frame) || frame <= 1 || frame > resolution) {
      frame <- ceiling(resolution / 2)
    }

    image(volume[, , frame], axes = FALSE, main = "1. Initial surface embed", asp = 1)

    image(volume_grew[, , frame] + volume[, , frame], axes =  FALSE, main = "2. Grow voxels", asp = 1)

    image(volume2[, , frame] + volume[, , frame], axes = FALSE,  main = "3. Bucket-fill+shrink", asp = 1)

  })

  # free up memory ~ 500MB
  # rm(volume, volume_shrunk, volume_filled, volume_grew)

  return(list(
    volume = volume2,
    IJK2RAS = IJK2RAS,
    resolution = resolution,
    fill = close_radius,
    close_radius = close_radius,
    close_radius_voxels = close_radius_voxels,
    pixdim = pixdim,
    inflate = inflate
  ))

}


