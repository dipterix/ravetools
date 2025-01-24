#' @title Sample '3D' volume in the world (anatomical \code{'RAS'}) space
#' @description
#' Low-level implementation to sample a '3D' volume into given orientation and
#' shape via a nearest-neighbor sampler.
#' @param x image (volume) to be sampled: \code{dim(x)} must have length of 3
#' @param new_dim target dimension, integers of length 3
#' @param vox2ras_old from volume index (column-row-slice) to \code{'RAS'}
#' (right-anterior-superior) transform: the volume index starts from 0 (C-style)
#' instead of 1 (R-style) to comply with \code{'NIfTI'} transform.
#' @param vox2ras_new the targeting transform from volume index to \code{'RAS'}
#' @param na_fill default numbers to fill if a pixel is out of bound; default is
#' \code{NA} or \code{as.raw(0)} if input \code{x} is raw type
#' @returns A newly sampled volume that aligns with \code{x} in the anatomical
#' \code{'RAS'} coordinate system. The underlying storage mode is the same as
#' \code{x}
#'
#' @examples
#'
#'
#' # up-sample and rotate image
#' x <- array(0, c(9, 9, 9))
#' x[4:6, 4:6, 4:6] <- 1
#' vox2ras <- matrix(nrow = 4, byrow = TRUE, c(
#'   0.7071, -0.7071, 0, 0,
#'   0.7071, 0.7071, 0, -5.5,
#'   0, 0, 1, -4,
#'   0, 0, 0, 1
#' ))
#'
#' new_vox2ras <- matrix(nrow = 4, byrow = TRUE, c(
#'   0, 0.5, 0, -4,
#'   0, 0, -0.5, 4,
#'   0.5, 0, 0, -4,
#'   0, 0, 0, 1
#' ))
#'
#' y <- resample_3d_volume(
#'   x,
#'   c(17, 17, 17),
#'   vox2ras_old = vox2ras,
#'   vox2ras_new = new_vox2ras,
#'   na_fill = 0
#' )
#'
#'
#' image(y[9,,])
#'
#'
#'
#' @export
resample_3d_volume <- function(x, new_dim, vox2ras_old, vox2ras_new = vox2ras_old, na_fill = NA) {
  verify_dim <- function(dm, name) {
    if(length(dm) < 3) {
      stop(sprintf("`resample_3d_volume`: %s must be a 3-dimensional volume.", name))
    }
    if(length(dm) > 3 && !all(dm[-c(1,2,3)] == 1)) {
      stop(
        sprintf("`resample_3d_volume`: %s has more than 3-dimensions and the extra ", name),
        "dimension has more than 1 component. Please make sure the volume dimension ",
        "has length of 3.")
    }
    dm[1:3]
  }
  verify_dim(dim(x), "Input `x`")
  new_dim <- verify_dim(new_dim, "Output")

  vox2ras_old <- as.matrix(vox2ras_old)
  if(!length(vox2ras_old) %in% c(12, 16) || ncol(vox2ras_old) != 4) {
    stop(
      "`resample_3d_volume`: `vox2ras_old` must be a 3x4 or 4x4 matrix ",
      "from the 0-indexing space to Right-Anterior-Superior coordinate system."
    )
  }

  vox2ras_new <- as.matrix(vox2ras_new)
  if(!length(vox2ras_new) %in% c(12, 16) || ncol(vox2ras_new) != 4) {
    stop(
      "`resample_3d_volume`: `vox2ras_new` must be a 3x4 or 4x4 matrix ",
      "from the 0-indexing space to Right-Anterior-Superior coordinate system."
    )
  }

  if(storage.mode(x) == "raw" ) {
    if(is.na(na_fill)) {
      na_fill <- as.raw(0)
    } else {
      storage.mode(na_fill) <- storage.mode(x)
    }
  } else {
    storage.mode(na_fill) <- storage.mode(x)
  }

  re <- resample3D(
    arrayDim = new_dim,
    fromArray = x,
    newVoxToWorldTransposed = t(vox2ras_new),
    oldVoxToWorldTransposed = t(vox2ras_old),
    na = na_fill
  )
  resampled <- re[[1]]

  vox2vox_orig <- re[[2]]
  dim(vox2vox_orig) <- c(4, 4)
  attr(resampled, "vox2ras") <- vox2ras_new
  attr(resampled, "vox2vox_orig") <- vox2vox_orig
  resampled
}
