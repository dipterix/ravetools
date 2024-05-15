#' @title Fill a volume cube based on water-tight surface
#' @author Zhengjia Wang
#' @description Create a cube volume (\code{256} 'voxels' on each margin), fill
#' in the 'voxels' that are inside of the surface.
#'
#' @param surface a surface mesh, can be mesh objects from \code{rgl} or
#' \code{freesurferformats} packages
#' @param inflate amount of 'voxels' to inflate on the final result; must be
#' a non-negative integer. A zero \code{inflate} value means the resulting
#' volume is tightly close to the surface
#' @param IJK2RAS volume 'IJK' (zero-indexed coordinate index) to
#' \code{'tkrRAS'} transform, default is automatically determined
#' leave it `NULL` if you don't know how to set it
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

ensure_mesh3d <- function(surface) {
  # check if surface is read from `freesurferformats`
  if (inherits(surface, "fs.surface")) {
    # surface <- freesurferformats::fs.surface.to.tmesh3d( surface )
    face <- t(surface$faces)
    if( min(face, na.rm = TRUE) == 0 ) {
      face <- face + 1
    }
    surface <- structure(list(
      vb = rbind(t(surface$vertices), 1),
      it = face), class = "mesh3d")
  }
  if (!inherits(surface, "mesh3d")) {
    # check if surface has names
    if(
      is.list(surface) &&
      "vb" %in% names(surface) &&
      is.matrix(surface$vb) &&
      nrow(surface$vb) %in% c(3L, 4L)
    ) {
      surface <- structure(surface, class = "mesh3d")
    } else {
      stop("ravetools:::ensure_mesh3d: `surface` must be a mesh3d object (package rgl) or a fs.surface object (package freesurferformats)")
    }
  }
  return(surface)
}

#' @rdname fill_surface
#' @export
fill_surface <- function(surface, inflate = 0, IJK2RAS = NULL, preview = FALSE,
                         preview_frame = 128) {
  # DIPSAUS DEBUG START
  # surface <- freesurferformats::read.fs.surface('~/Dropbox (PENN Neurotrauma)/RAVE/Samples/raw/PAV010/rave-imaging/fs/surf/rh.pial')
  # IJK2RAS <- NULL
  # verbose <- TRUE
  # preview <- TRUE
  # output_format <- "rgl"
  # inflate <- 3

  inflate <- as.integer(inflate)
  if(inflate < 0) { inflate <- 0 }

  # setup
  preview2D <- function(expr) {
    expr <- substitute(expr)
    if(interactive() && preview) {
      eval(expr)
    }
  }

  # Initialize the parameters
  resolution <- 256L
  if(length(IJK2RAS) != 16L) {
    IJK2RAS <- matrix(
      nrow = 4, byrow = TRUE,
      c(-1, 0, 0, resolution/2,
        0, 0, 1, -resolution/2,
        0, -1, 0, resolution/2,
        0, 0, 0, 1))
  }

  # Make sure the surface is a `mesh3d` object
  surface_orig <- ensure_mesh3d(surface)

  # Transform surface to IJK voxel space so can be fitted into a volume
  # Remember IJK starts from 0
  surface <- surface_orig
  surface$vb <- solve(IJK2RAS) %*% surface$vb

  # Get maximum voxel numbers that can be used to fill in
  vertex_range <- apply(surface$vb, 1, range)
  max_fill <- floor(min(c(
    vertex_range[1, 1:3],
    resolution - vertex_range[2, 1:3],
    17
  )) - 2L)

  # Creating a volume
  volume <- array(0L, dim = rep(resolution, 3))

  # embed the surface in volume space
  surface_index <- round(surface$vb[c(1,2,3), ])
  surface_index <- surface_index[1, ] + surface_index[2, ] * resolution +
    surface_index[3, ] * (resolution^2) + 1
  volume[surface_index] <- 1L

  # Grow the volume by ~15mm and shrink back (~12mm).
  # This step connects the
  # segmented voxels into a shell that is water-tight
  volume_grew <- grow_volume(volume, max_fill)

  # bucket-fill from the corner, bucketFillVolume is in-place
  volume_filled <- bucketFillVolume(volume_grew + 0L, x = 1, y = 1, z = 1, fill = 1L)

  # Fill the voxels within the surface and shrink
  shrink_amount <- max_fill - inflate
  volume2 <- 1L - grow_volume(volume_filled - volume_grew, shrink_amount)

  # preview
  preview2D({
    oldPar <- graphics::par(c("mfrow", "mar"))
    graphics::par(mfrow = c(1,3), mar = c(3.1, 0.1, 3.1, 0.1))
    on.exit({ do.call(graphics::par, oldPar) })

    frame <- as.integer(preview_frame)
    if(is.na(frame) || frame <= 1 || frame > resolution) {
      frame <- ceiling(resolution / 2)
    }

    image(volume[,,frame], axes = FALSE, main = "1. Initial surface embed", asp = 1)

    image(volume_grew[,,frame] + volume[,,frame], axes = FALSE, main = "2. Grow voxels", asp = 1)

    image(volume2[,,frame] + volume[,,frame], axes = FALSE, main = "3. Bucket-fill+shrink", asp = 1)

  })

  # free up memory ~ 500MB
  # rm(volume, volume_shrunk, volume_filled, volume_grew)

  return(list(
    volume = volume2,
    IJK2RAS = IJK2RAS,
    fill = max_fill,
    inflate = inflate
  ))

}


