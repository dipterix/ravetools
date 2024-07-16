meshintegrity <- function(mesh, facecheck = FALSE, normcheck = FALSE) {

  mesh <- ensure_mesh3d(mesh)

  ##check vertices
  if (!is.null(mesh$vb) && is.matrix(mesh$vb)) {
    vdim <- dim(mesh$vb)
    if (vdim[2] < 1) {
      stop("mesh has no vertices")
    }
    if (!vdim[1] %in% c(3, 4)) {
      stop("vertices have invalid dimensionality")
    }
    if(!identical(storage.mode(mesh$vb), "double")) {
      storage.mode(mesh$vb) <- "double"
    }
    if ( anyNA(mesh$vb) ) {
      stop("vertex coords need to be numeric values")
    }
  } else {
    stop("mesh has no/invalid vertices")
  }
  if (is.matrix(mesh$it)){
    if (ncol(mesh$it) == 0)
      mesh$it <- NULL
  }
  if (!is.null(mesh$it) && is.matrix(mesh$it)) {
    itdim <- dim(mesh$it)
    if(itdim[1] != 3) {
      stop("only triangular faces are valid")
    }

    if(!identical(storage.mode(mesh$it), "integer")) {
      storage.mode(mesh$it) <- "integer"
    }
    if( anyNA(mesh$it) ) {
      stop("face indices need to be integer values")
    }

    itrange <- range(mesh$it)
    if (itrange[1] < 1 || itrange[2] > vdim[2]) {
      stop("faces reference non-existent vertices")
    }
  } else if (facecheck) {
    stop("mesh has no triangular face indices")
  }

  if (normcheck) {
    if (!is.null(mesh$normals) && is.matrix(mesh$normals)) {
      ndim <- dim(mesh$normals)
      if(!prod(ndim == vdim)) {
        stop("normals must be of same dimensionality as vertices")
      }

      if ( anyNA(mesh$normals) ) {
        stop("normal coords need to be numeric values")
      }
    } else {
      stop("mesh has no vertex normals")
    }
  }
  return(mesh)
}

cSizeMesh <- function(mesh) {
  x <- t(mesh$vb[1:3, ])
  X <- scale(x, scale = FALSE)
  y <- sqrt(sum(as.vector(X) ^ 2))
  return(y)
}

pcax <- function(mesh) {
  x <- t(mesh$vb[1:3, ])
  pc <- 3 * stats::prcomp(x, retx = FALSE)$sdev[1]
  return(pc)
}

meshOff <- function(x, offset) {
  x <- vcg_update_normals(x)
  x$vb[1:3, ] <- x$vb[1:3, ] + offset * x$normals[1:3, ]
  return(x)
}

bbox <- function(x) {
  # Bounding box

  bbox <- apply(x$vb[1:3, , drop = FALSE], 1L, range)
  bbox <- expand.grid(bbox[, 1], bbox[, 2], bbox[, 3])
  dia <- max(stats::dist(bbox))
  return(list(bbox = bbox, diag = dia))
}

checkFaceOrientation <- function(x, offset = NULL) {
  if (is.null(offset)) {
    offset <- pcax(x) / 15
  }
  out <- TRUE
  xoff <- meshOff(x, offset)
  cx <- cSizeMesh(x)
  cxoff <- cSizeMesh(xoff)
  if (cx > cxoff) {
    out <- FALSE
  }
  return(out)
}

invertFaces <- function (mesh) {
  mesh$it <- mesh$it[c(3, 2, 1), ]
  mesh <- vcg_update_normals(mesh)
  return(mesh)
}

#' @title Update vertex normal
#' @param mesh triangular mesh or a point-cloud (matrix of 3 columns)
#' @param weight method to compute per-vertex normal vectors: \code{"area"}
#' weighted average of surrounding face normal, or \code{"angle"} weighted
#' vertex normal vectors.
#' @param pointcloud integer vector of length 2: containing optional
#' parameters for normal calculation of point clouds; the first entry
#' specifies the number of neighboring points to consider; the second
#' entry specifies the amount of smoothing iterations to be performed.
#' @param verbose whether to verbose the progress
#' @returns A \code{'mesh3d'} object with normal vectors.
#' @examples
#'
#' if(is_not_cran()) {
#'
#' # Prepare mesh with no normal
#' data("left_hippocampus_mask")
#' mesh <- vcg_isosurface(left_hippocampus_mask)
#' mesh$normals <- NULL
#'
#' # Start: examples
#' new_mesh <- vcg_update_normals(mesh, weight = "angle",
#'                                pointcloud = c(10, 10))
#'
#' rgl_view({
#'   rgl_call("mfrow3d", 1, 2)
#'   rgl_call("shade3d", mesh, col = 2)
#'
#'   rgl_call("next3d")
#'   rgl_call("shade3d", new_mesh, col = 2)
#' })
#' }
#'
#'
#' @export
vcg_update_normals <- function(
    mesh,
    weight = c("area", "angle"),
    pointcloud = c(10, 0),
    verbose = FALSE
) {
  weight <- match.arg(weight)
  if (length(pointcloud) != 2) {
    stop("pointcloud must be an integer vector of length 2")
  }

  if( weight == "area" ) {
    type <- 0L
  } else {
    type <- 1L
  }
  if ( is.matrix(mesh)) {
    tmp <- list()
    tmp$vb <- rbind(t(mesh), 1)
    mesh <- tmp
    class(mesh) <- "mesh3d"
  }
  mesh <- meshintegrity(mesh)
  vb <- mesh$vb[1:3, , drop = FALSE]
  if (!is.matrix(vb)) {
    stop("mesh has no vertices")
  }
  it <- mesh$it - 1L
  normals <- vcgUpdateNormals(vb, it, type, pointcloud, !verbose)
  mesh$normals <- rbind(normals, 1)
  return(mesh)
}



#' @name vcg_smooth
#' @title Implicitly smooth a triangular mesh
#' @description
#' Applies smoothing algorithms on a triangular mesh.
#'
#' @param mesh triangular mesh stored as object of class 'mesh3d'.
#' @param use_mass_matrix logical: whether to use mass matrix to keep the mesh
#' close to its original position (weighted per area distributed on vertices);
#' default is \code{TRUE}
#' @param fix_border logical: whether to fix the border vertices of the mesh;
#' default is \code{FALSE}
#' @param use_cot_weight logical: whether to use cotangent weight; default is
#' \code{FALSE} (using uniform 'Laplacian')
#' @param laplacian_weight numeric: weight when \code{use_cot_weight} is \code{FALSE};
#' default is \code{1.0}
#' @param degree integer: degrees of 'Laplacian'; default is \code{1}
#' @param type method name of explicit smooth, choices are \code{'taubin'},
#' \code{'laplace'}, \code{'HClaplace'}, \code{'fujiLaplace'},
#' \code{'angWeight'}, \code{'surfPreserveLaplace'}.
#' @param iteration number of iterations
#' @param lambda In \code{vcg_smooth_implicit}, the amount of smoothness,
#' useful only if \code{use_mass_matrix} is \code{TRUE}; default is \code{0.2}.
#' In \code{vcg_smooth_explicit}, parameter for \code{'taubin'} smoothing.
#' @param mu parameter for \code{'taubin'} explicit smoothing.
#' @param delta parameter for scale-dependent 'Laplacian' smoothing or
#' maximum allowed angle (in 'Radian') for deviation between surface preserving
#' 'Laplacian'.
#' @returns An object of class "mesh3d" with:
#' \item{\code{vb}}{vertex coordinates}
#' \item{\code{normals}}{vertex normal vectors}
#' \item{\code{it}}{triangular face index}
#' @examples
#'
#' if(is_not_cran()) {
#'
#' # Prepare mesh with no normals
#' data("left_hippocampus_mask")
#'
#' # Grow 2mm on each direction to fill holes
#' volume <- grow_volume(left_hippocampus_mask, 2)
#'
#' # Initial mesh
#' mesh <- vcg_isosurface(volume)
#'
#' # Start: examples
#' rgl_view({
#'   rgl_call("mfrow3d", 2, 4)
#'   rgl_call("title3d", "Naive ISOSurface")
#'   rgl_call("shade3d", mesh, col = 2)
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "Implicit Smooth")
#'   rgl_call("shade3d", col = 2,
#'            x = vcg_smooth_implicit(mesh, degree = 2))
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "Explicit Smooth - taubin")
#'   rgl_call("shade3d", col = 2,
#'            x = vcg_smooth_explicit(mesh, "taubin"))
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "Explicit Smooth - laplace")
#'   rgl_call("shade3d", col = 2,
#'            x = vcg_smooth_explicit(mesh, "laplace"))
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "Explicit Smooth - angWeight")
#'   rgl_call("shade3d", col = 2,
#'            x = vcg_smooth_explicit(mesh, "angWeight"))
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "Explicit Smooth - HClaplace")
#'   rgl_call("shade3d", col = 2,
#'            x = vcg_smooth_explicit(mesh, "HClaplace"))
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "Explicit Smooth - fujiLaplace")
#'   rgl_call("shade3d", col = 2,
#'            x = vcg_smooth_explicit(mesh, "fujiLaplace"))
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "Explicit Smooth - surfPreserveLaplace")
#'   rgl_call("shade3d", col = 2,
#'            x = vcg_smooth_explicit(mesh, "surfPreserveLaplace"))
#' })
#'
#' }
#'
#' @export
vcg_smooth_implicit <- function(
    mesh, lambda = 0.2, use_mass_matrix = TRUE, fix_border = FALSE,
    use_cot_weight = FALSE, degree = 1L, laplacian_weight = 1.0
) {
  mesh <- meshintegrity(mesh)
  smooth_quality <- FALSE
  vb <- mesh$vb[1:3, , drop = FALSE]
  it <- (mesh$it - 1L)
  storage.mode(it) <- "integer"
  lambda <- as.double(lambda)[[1]]
  laplacian_weight <- as.double(laplacian_weight)[[1]]
  use_mass_matrix <- as.logical(use_mass_matrix)[[1]]
  fix_border <- as.logical(fix_border)[[1]]
  use_cot_weight <- as.logical(use_cot_weight)[[1]]
  smooth_quality <- as.logical(smooth_quality)[[1]]
  degree <- as.integer(degree)[[1]]
  tmp <- vcgSmoothImplicit(vb, it, lambda, use_mass_matrix, fix_border, use_cot_weight, degree, laplacian_weight, smooth_quality)
  mesh$vb[1:3, ] <- tmp$vb
  mesh$normals <- rbind(tmp$normals, 1)
  mesh$it <- tmp$it
  mesh
}

#' @rdname vcg_smooth
#' @export
vcg_smooth_explicit <- function(
  mesh, type = c("taubin", "laplace", "HClaplace", "fujiLaplace", "angWeight", "surfPreserveLaplace"),
  iteration = 10, lambda = 0.5, mu = -0.53, delta = 0.1
) {
  mesh <- meshintegrity(mesh)
  type <- match.arg(type)
  type <- substring(type[1],1L,1L)
  vb <- mesh$vb[1:3, , drop = FALSE]
  it <- (mesh$it - 1L)
  storage.mode(it) <- "integer"
  method <- 0
  if (type == "l" || type == "L") {
    method <- 1
  } else if (type == "H" || type == "h") {
    method <- 2
  } else if (type == "f" || type == "F") {
    method <- 3
  } else if (type == "a" || type == "A") {
    method <- 4
  } else if (type == "s" || type == "S") {
    method <- 5
  }

  stopifnot(is.integer(it))

  tmp <- vcgSmooth(vb, it, iteration, method, lambda, mu, delta)
  mesh$vb[1:3,] <- tmp$vb
  mesh$normals <- rbind(tmp$normals, 1)
  mesh$it <- tmp$it
  invisible(meshintegrity(mesh))
}

#' @title Compute volume for manifold meshes
#' @param mesh triangular mesh of class \code{'mesh3d'}
#' @returns The numeric volume of the mesh
#'
#' @examples
#'
#' # Initial mesh
#' mesh <- vcg_sphere()
#'
#' vcg_mesh_volume(mesh)
#'
#' @export
vcg_mesh_volume <- function(mesh) {
  vcgVolume( meshintegrity(mesh) )
}

#' Simple 3-dimensional sphere mesh
#' @param sub_division density of vertex in the resulting mesh
#' @param normals whether the normal vectors should be calculated
#' @returns A \code{'mesh3d'} object
#' @examples
#'
#' vcg_sphere()
#'
#' @export
vcg_sphere <- function(sub_division = 3L, normals = TRUE) {
  vcgSphere(sub_division, normals)
}

#' @title Create surface mesh from 3D-array
#' @description
#' Create surface from 3D-array using marching cubes algorithm
#'
#' @param volume a volume or a mask volume
#' @param threshold_lb lower-bound threshold for creating the surface; default
#' is \code{0}
#' @param threshold_ub upper-bound threshold for creating the surface; default
#' is \code{NA} (no upper-bound)
#' @param vox_to_ras a \code{4x4} \code{'affine'} transform matrix indicating the
#' 'voxel'-to-world transform.
#'
#' @returns A triangular mesh of class \code{'mesh3d'}
#'
#' @examples
#'
#'
#' if(is_not_cran()) {
#'
#' library(ravetools)
#' data("left_hippocampus_mask")
#'
#' mesh <- vcg_isosurface(left_hippocampus_mask)
#'
#'
#' rgl_view({
#'
#'   rgl_call("mfrow3d", 1, 2)
#'
#'   rgl_call("title3d", "Direct ISOSurface")
#'   rgl_call("shade3d", mesh, col = 2)
#'
#'   rgl_call("next3d")
#'   rgl_call("title3d", "ISOSurface + Implicit Smooth")
#'
#'   rgl_call("shade3d",
#'            vcg_smooth_implicit(mesh, degree = 2),
#'            col = 3)
#' })
#'
#' }
#' @export
vcg_isosurface <- function(
    volume,
    threshold_lb = 0, threshold_ub = NA,
    vox_to_ras = diag(c(-1, -1, 1, 1))
) {
  if (!length(volume) || length(dim(volume)) != 3) {
    stop("vcg_isosurface: 3D non-empty `volume` array needed")
  }

  if(is.na(threshold_lb)) { threshold_lb <- 0 }
  sel <- volume > threshold_lb
  dimnames(sel) <- NULL

  if(!is.na(threshold_ub)) {
    sel <- sel & volume < threshold_ub
  }

  mesh <- structure(vcgIsoSurface( sel, 0.5 ), class = "mesh3d")

  # voxel (0-indexed) to RAS
  mesh$vb <- (vox_to_ras %*% rbind(mesh$vb, 1))[seq_len(3), ]

  if (!checkFaceOrientation(mesh)) {
    mesh <- invertFaces(mesh)
  }
  return(mesh)
}



#' Sample a surface mesh uniformly
#' @param x surface
#' @param voxel_size 'voxel' size for space 'discretization'
#' @param offset offset position shift of the new surface from the input
#' @param discretize whether to use step function (\code{TRUE}) instead of
#' linear interpolation (\code{FALSE}) to calculate the position of the
#' intersected edge of the marching cube; default is \code{FALSE}
#' @param multi_sample whether to calculate multiple samples for more accurate
#' results (at the expense of more computing time) to remove artifacts; default
#' is \code{FALSE}
#' @param absolute_distance whether an unsigned distance field should be
#' computed. When set to \code{TRUE}, non-zero offsets is to be set, and
#' double-surfaces will be built around the original surface, like a sandwich.
#' @param merge_clost whether to merge close vertices; default is \code{TRUE}
#' @param verbose whether to verbose the progress; default is \code{TRUE}
#' @returns A triangular mesh of class \code{'mesh3d'}
#' @examples
#'
#' sphere <- vcg_sphere()
#' mesh <- vcg_uniform_remesh(sphere, voxel_size = 0.45)
#'
#' if(is_not_cran()) {
#'
#' rgl_view({
#'
#'   rgl_call("mfrow3d", 1, 2)
#'
#'   rgl_call("title3d", "Input")
#'   rgl_call("wire3d", sphere, col = 2)
#'   rgl_call("next3d")
#'
#'   rgl_call("title3d", "Re-meshed to 0.1mm edge distance")
#'   rgl_call("wire3d", mesh, col = 3)
#' })
#'
#' }
#'
#' @export
vcg_uniform_remesh <- function(
    x, voxel_size = NULL, offset = 0, discretize = FALSE,
    multi_sample = FALSE, absolute_distance = FALSE, merge_clost = FALSE, verbose = TRUE) {

  x <- meshintegrity(mesh = x, facecheck = TRUE)
  if (is.null( voxel_size )) {
    voxel_size <- bbox(x)$dia / 50
  }
  vb <- x$vb
  it <- x$it - 1L
  out <- structure(
    vcgUniformResample( vb, it, voxel_size, offset, discretize,
                        multi_sample, absolute_distance, merge_clost, !verbose ),
    class = "mesh3d"
  )
  return(meshintegrity(out))
}
