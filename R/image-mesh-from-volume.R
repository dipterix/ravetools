mesh_volume <- function(mesh) {
  suppressWarnings({
    tryCatch({
      return(vcg_mesh_volume(mesh))
    }, error = function(e) {
      NA_real_
    })
  })
}

#' @title Generate 3D mesh surface from volume data
#' @description This function is soft-deprecated. Please use
#' \code{\link{vcg_mesh_volume}}, \code{\link{vcg_uniform_remesh}}, and
#' \code{\link{vcg_smooth_explicit}} or \code{\link{vcg_smooth_implicit}}.
#' @param volume 3-dimensional volume array
#' @param output_format resulting data format, choices are \code{'rgl'} and
#' \code{'freesurfer'}
#' @param IJK2RAS volume 'IJK' (zero-indexed coordinate index) to
#' \code{'tkrRAS'} transform, default is automatically determined
#' @param threshold threshold used to create volume mask; the surface will be
#' created to fit the mask boundaries
#' @param verbose whether to verbose the progress
#' @param remesh whether to re-sample the mesh using \code{\link{vcg_uniform_remesh}}
#' @param remesh_voxel_size,remesh_multisample,remesh_automerge see
#' arguments in \code{\link{vcg_uniform_remesh}}
#' @param smooth whether to smooth the mesh via \code{\link{vcg_smooth_explicit}}
#' @param smooth_lambda,smooth_delta,smooth_method see \code{\link{vcg_smooth_explicit}}
#' @returns A \code{'mesh3d'} surface if \code{output_format} is 'rgl', or
#' \code{'fs.surface'} surface otherwise.
#' @examples
#'
#'
#' volume <- array(0, dim = c(8,8,8))
#' volume[4:5, 4:5, 4:5] <- 1
#'
#' graphics::image(x = volume[4,,])
#'
#' # you can use rgl::wire3d(mesh) to visualize the mesh
#' mesh <- mesh_from_volume(volume, verbose = FALSE)
#'
#'
#' @export
mesh_from_volume <- function(
    volume, output_format = c("rgl", "freesurfer"),
    IJK2RAS = NULL, threshold = 0, verbose = TRUE,
    remesh = TRUE, remesh_voxel_size = 1,
    remesh_multisample = TRUE, remesh_automerge = TRUE,
    smooth = FALSE, smooth_lambda = 10, smooth_delta = 20,
    smooth_method = "surfPreserveLaplace"
) {

  output_format <- match.arg(output_format)

  debug <- function(..., appendLF = TRUE) {
    if (verbose) {
      appendLF <- ifelse(appendLF, "\n", "")
      cat(..., appendLF)
    }
  }

  volume <- as.array(volume)
  dm <- dim(volume)[c(1,2,3)]
  dim(volume) <- dm

  if(length(IJK2RAS) != 16L) {
    IJK2RAS <- matrix(
      nrow = 4, byrow = TRUE,
      c(-1, 0, 0, dm[[1]]/2,
        0, 0, 1, -dm[[3]]/2,
        0, -1, 0, dm[[2]]/2,
        0, 0, 0, 1))
  }

  mesh <- vcg_isosurface(volume, threshold_lb = threshold, vox_to_ras = IJK2RAS)

  debug(sprintf("The initial reconstructed surface volume is %.1f mm^3", mesh_volume(mesh)))

  if( remesh ) {
    mesh <- vcg_uniform_remesh(
      mesh, voxel_size = remesh_voxel_size, multi_sample = remesh_multisample,
      merge_clost = remesh_automerge, verbose = verbose)
    debug(sprintf("The re-meshed surface volume is %.1f mm^3", mesh_volume(mesh)))
  }
  if( smooth ) {
    mesh <- vcg_smooth_explicit(
      mesh = mesh, type = smooth_method,
      lambda = smooth_lambda, delta = smooth_delta
    )
    debug(sprintf("The smoothed surface volume is %.1f mm^3", mesh_volume(mesh)))
  }

  switch(
    output_format,
    "freesurfer" = {
      face_index <- t(mesh$it)
      face_index_start <- min(face_index)
      face_index <- face_index - (face_index_start - 1L)

      vertices <- t(mesh$vb[c(1, 2, 3), ])
      mesh <- structure(
        list(
          mesh_face_type = "tris",
          vertices = vertices,
          faces = face_index
        ),
        class = c("fs.surface", "list")
      )
    },{}
  )
  mesh

}
