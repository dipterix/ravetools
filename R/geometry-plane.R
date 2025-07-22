#' @title Create a two-dimensional plane in three dimensional space
#' @param width,height width and height of the plane, must not be \code{NA}
#' @param shape length of two to indicate the number of vertices along width
#' and height, default is only \code{c(2, 2)} (2 vertices each side, hence
#' one grid)
#' @returns A triangular mesh of class \code{'mesh3d'}
#' @examples
#'
#' plane <- plane_geometry(5, 10, c(12, 22))
#'
#' if(FALSE) {
#'
#'   rgl_view({
#'
#'     rgl_call("shade3d", plane, col = 3)
#'     rgl_call("wire3d", plane, col = 1)
#'
#'   })
#'
#' }
#'
#' @export
plane_geometry <- function(width = 1, height = 1, shape = c(2, 2)) {

  stopifnot2(isTRUE(is.finite(width)), msg = "plane_geometry: `width` must be a valid number")
  stopifnot2(isTRUE(is.finite(height)), msg = "plane_geometry: `height` must be a valid number")

  width_half <- width / 2
  height_half <- height / 2

  grid_size <- shape - 1

  grid_size <- floor(grid_size)
  if(length(grid_size) == 1) {
    grid_size <- c(grid_size, grid_size)
  } else {
    grid_size <- grid_size[c(1, 2)]
  }

  grid_x1 <- grid_size[[1]] + 1
  grid_y1 <- grid_size[[2]] + 1

  segment_size <- c(width, height) / grid_size

  ix <- seq(0, grid_size[[1]])
  vertices <- sapply(seq(0, grid_size[[2]]), function(iy) {
    y <- iy * segment_size[[2]] - height_half

    x <- ix * segment_size[[1]] - width_half

    rbind(x, -y, 0)
  }, simplify = TRUE, USE.NAMES = FALSE)

  vertices <- matrix(as.vector(vertices), nrow = 3L, byrow = FALSE)

  ix <- seq_len(grid_size[[1]])
  indices <- sapply(seq_len(grid_size[[2]]), function(iy) {
    a <- ix + grid_x1 * (iy - 1)
    b <- ix + grid_x1 * iy
    c <- ix + 1 + grid_x1 * iy
    d <- ix + 1 + grid_x1 * (iy - 1)
    rbind(a, b, d, b, c, d)
  })
  indices <- matrix(as.vector(indices), nrow = 3L, byrow = FALSE)

  # return rgl-style surface
  structure(list(vb = vertices, it = indices), class = "mesh3d")
}

#' @title Project plane to a surface
#' @description
#' Project a two-dimensional plane (such as \code{'ECoG'} grid) to a
#' three-dimensional surface while preserving the order
#' @param target target surface to be projected to, must be object that can
#' be converted to \code{'mesh3d'} (\code{'rgl'} surface), for example,
#' \code{'fs.surface'} (from \code{'freesurferformat'} package) or
#' \code{'ieegio_surface'} from \code{'ieegio'} package.
#' @param width,height width and height of the plane in world space (for
#' \code{'ECoG'} grid, the unit is millimeter)
#' @param shape vector of two integers: the first element is the number of
#' vertices (or electrode contacts) along \code{'width'} direction; the second
#' element is the number of vertices along \code{'height'} direction. The total
#' number of vertices of the plane will be \code{prod(shape)}. Notice
#' @param initial_positions a \code{shape[[1]] x shape[[2]] x 3} array or a
#' \code{n x 3} matrix, where \code{n} is \code{prod(shape)}, the number of
#' vertices indicating the initial vertex positions of the plane
#' @param translate_first whether to translate the plane first if the
#' plane center is far from the surface; default is \code{FALSE}; set to
#' \code{TRUE} for a warm start
#' @param diagnostic whether to plot diagnostic figures showing the morphing
#' progress.
#' @returns The projected vertex locations, same order as \code{initial_positions}.
#'
#' @examples
#'
#'
#' # Construct target surface
#'
#' sphere <- vcg_sphere()
#' target <- structure(
#'   class = "mesh3d",
#'   list(
#'     vb = cbind(
#'       sphere$vb[1:3, ] - c(0.8, 0, 0),
#'       sphere$vb[1:3, ] + c(0.8, 0, 0)
#'     ),
#'     it = cbind(
#'       sphere$it[1:3, ],
#'       sphere$it[1:3, ] + ncol(sphere$vb)
#'     )
#'   )
#' )
#' n_surfverts <- ncol(target$vb)
#'
#' plane <- plane_geometry(width = 3, height = 3, shape = c(30, 30))
#' plane$vb <- plane$vb[1:3, , drop = FALSE] + c(0, 0, 2)
#' n_contacts <- ncol(plane$vb)
#'
#' # First plot
#' x <- t(cbind(target$vb, plane$vb))
#' colnames(x) <- c('x', 'y', 'z')
#' graphics::pairs(
#'   x = x, asp = 1,
#'   col = c(
#'     rep("black", n_surfverts),
#'     rep("green", n_contacts)
#'   ),
#'   pch = c(
#'     rep(46, n_surfverts),
#'     rep(20, n_contacts)
#'   )
#' )
#'
#' projected <- project_plane(
#'   target = target, width = 3, height = 3, shape = c(30, 30),
#'   initial_positions = t(plane$vb),
#'   translate_first = TRUE, diagnostic = FALSE
#' )
#'
#' y <- rbind(x, projected)
#' graphics::pairs(
#'   x = y, asp = 1,
#'   col = c(
#'     rep("black", ncol(target$vb)),
#'     rep("green", n_contacts),
#'     rep("red", n_contacts)
#'   ),
#'   pch = c(
#'     rep(46, n_surfverts),
#'     rep(1, n_contacts),
#'     rep(20, n_contacts)
#'   )
#' )
#'
#'
#' @export
project_plane <- function(
    target, width, height, shape, initial_positions,
    translate_first = TRUE, diagnostic = FALSE) {

  target <- ensure_mesh3d(target)
  plane <- plane_geometry(width = width, height = height, shape = shape)

  # set initial electrode locations
  n_contacts <- ncol(plane$vb)
  if(n_contacts == 1) {
    initial_positions <- matrix(initial_positions, nrow = 1, ncol = 3)
  } else {
    if(!is.matrix(initial_positions) && is.array(initial_positions)) {
      dm <- dim(initial_positions)
      dim(initial_positions) <- c(dm[[1]] * dm[[2]], dm[[3]])
    }
    stopifnot2(ncol(initial_positions) == 3, msg = "project_plane: `initial_positions` must be a matrix with 3 columns")
    stopifnot2(nrow(initial_positions) == n_contacts, msg = "project_plane: number of `initial_positions` (row count) must equal to `shape[1] x shape[2]`")
    initial_positions <- as.matrix(initial_positions)
  }

  # keep a copy of model positions
  plane$vb[1:3, ] <- t(initial_positions)
  plane_orig <- plane
  # update normals
  plane <- vcg_update_normals(plane)

  # register utility functions
  diagnose_plot <- function() {
    if(!diagnostic) { return() }

    x <- t(cbind(plane_orig$vb[1:3, , drop = FALSE], plane$vb[1:3, , drop = FALSE]))
    colnames(x) <- c("x", "y", "z")
    col <- c(rep("green", n_contacts), rep("red", n_contacts))
    pch <- c(rep(1, n_contacts), rep(20, n_contacts))

    graphics::pairs(x, col = col, pch = pch,
                    asp = 1,
                    main = sprintf("Diagnostic plot (green=init,red=current)\ndim: x=%.1f y=%.1f, shape: %dx%d",
                                   width, height, shape[[1]], shape[[2]]))
  }

  nearest_nodes <- function(v_ind, radius) {
    dist <- sqrt(colSums((plane_orig$vb - plane_orig$vb[, v_ind])^2))
    which(dist <= radius)
  }

  # first projection - project the whole grid center to the surface
  if( translate_first ) {

    # compute the distances from contacts to the surface
    search_results <- vcg_kdtree_nearest(target, query = plane, k = 1)
    if(mean(search_results$distance) > 1) {

      # translate is needed when all contacts are atleast 1mm away from the surface

      average_normals <- rowMeans(plane$normals[1:3, , drop = FALSE])
      # project the plane along
      raycaster <- vcg_raycaster(
        target,
        ray_origin = plane$vb[1:3, , drop = FALSE],
        ray_direction = average_normals,
        both_sides = TRUE
      )

      # Find within contacts with ray-casting intersection
      idx_has_intersection <- which(raycaster$has_intersection)

      if(length(idx_has_intersection)) {

        contact_idx <- floor((shape[[2]] - 1) / 2) * shape[[1]] + round(shape[[1]] / 2)
        if(!contact_idx %in% idx_has_intersection) {
          idx2 <- which.min(abs(raycaster$distance)[idx_has_intersection])
          contact_idx <- idx_has_intersection[idx2][[1]]
        }
        dist <- raycaster$distance[contact_idx]
        proj_dir <- raycaster$ray_direction[, contact_idx, drop = TRUE] * dist

        # move the plane along the direction
        plane$vb[1:3, ] <- plane$vb[1:3, ] + proj_dir

      }

    }

    diagnose_plot()

  }

  # Total number of iterations are n_contacts x 5
  # 1, 2. Hebbian learning (random contacts)
  # 3, 4. iterate all contacts with large radius
  # 5. refine with small radius

  index_sampler <- c(sample(n_contacts), sample(n_contacts), sample(n_contacts), sample(n_contacts))
  lambda_fun <- function(iter) {
    # 0.2, 0.1, 0.2, 0.2, 0.2
    if(iter <= n_contacts) { return(0.2) }
    if(iter <= 2 * n_contacts) { return(0.1) }
    if(iter <= 3 * n_contacts) { return(0.2) }
    if(iter <= 4 * n_contacts) { return(0.2) }
    if(iter <= 5 * n_contacts) { return(0.2) }
  }

  plane_diag <- sqrt(sum(c(width, height)^2))
  contact_diag <- sqrt(sum(c(width / (shape[[1]] + 1), height / (shape[[1]] + 1))^2)) + 1e-5
  radius_ratio <- function(iter) {
    if(iter <= min(n_contacts, 100)) { return( 0.5 ) }
    if(iter <= n_contacts) { return(0.3) }
    if(iter <= 2 * n_contacts) { return(0.2) }
    if(iter <= 3 * n_contacts) { return(0.2) }
    if(iter <= 4 * n_contacts) { return(0.3) }
    if(iter <= 5 * n_contacts) { return(0.1) }
  }
  radius_fun <- function(iter) {
    max(radius_ratio(iter) * plane_diag, contact_diag)
  }

  # iterate
  for(iter in seq_len(n_contacts * 5)) {

    lambda <- lambda_fun(iter)
    radius <- radius_fun(iter)

    # for every grand iteration, update the projection twice
    if(iter %% floor(n_contacts / 2) == 1) {
      search_results <- vcg_kdtree_nearest(target, query = plane, k = 1)
    }

    # project direction (to nearest location on mesh)
    proj_dir <- target$vb[1:3, search_results$index, drop = FALSE] - plane$vb[1:3, , drop = FALSE]
    proj_dist <- sqrt(colSums(proj_dir^2))
    proj_dist[proj_dist <= 0] <- 1
    proj_dir <- t(t(proj_dir) / proj_dist)

    # dot product to electrode normal - we want to project electrodes,
    # but we do not want to bend the electrode too much
    dot_prod <- abs(colSums(proj_dir * plane$normals[1:3, ]))

    if(iter <= n_contacts) {
      # These are easy ones: the contacts are close to the mesh
      # or the projection aligns with the normals
      if(stats::runif(1) > 0.85) {
        # have prob to project distanced contacts
        min_idx <- sample(
          c(which(search_results$distance >= stats::quantile(search_results$distance, 0.9))),
          size = 1
        )
      } else {
        min_idx <- sample(
          c(which(dot_prod > 0.95), which.min(search_results$distance)),
          size = 1
        )
      }
    } else {
      min_idx <- index_sampler[iter - n_contacts]
    }

    neib <- nearest_nodes(v_ind = min_idx, radius = radius)
    surface_points <- target$vb[1:3, search_results$index[neib], drop = FALSE]
    update_vec <- rowMeans(surface_points - plane$vb[1:3, neib, drop = FALSE])

    if(dot_prod[min_idx] <= 0.86) {
      # the projection direction and plane normals have at least 30-degree dis-agreement
      # use the plane normal instead
      update_vec2 <- plane$normals[1:3, min_idx] * sqrt(sum(update_vec^2))
      if(sum(update_vec2 * update_vec) <= 0) {
        update_vec2 <- -update_vec2
      }
      update_vec <- update_vec2
    }

    # update_nodes <- nearest_nodes(min_idx, radius = radius)
    update_nodes <- neib

    # update plane
    plane$vb[1:3, update_nodes] <- plane$vb[1:3, update_nodes] + update_vec * lambda

    # smooth
    plane <- vcg_update_normals(plane, weight = "area")
    if(iter %% n_contacts == 0) {
      plane <- vcg_smooth_explicit(plane, type = "taubin", iteration = 2)
      plane <- vcg_update_normals(plane, weight = "area")
    }

    if( diagnostic && iter %% 100 == 0) {
      qdist <- stats::quantile(proj_dist, c(0.5, 1))
      dmsg <- sprintf("\r%d, current_dist=%.2f, media_dist=%.2f, max_dist=%.2f",
                      iter, proj_dist[min_idx], qdist[[1]], qdist[[2]])
      diagnose_plot()
      message(dmsg, appendLF = FALSE)
    }

  }

  # Finally, ray-caster
  # determine the hemisphere
  plane_center <- rowMeans(plane$vb)[1:3]
  target_center <- rowMeans(target$vb)[1:3]
  dir <- plane_center - target_center
  vec3 <- new_vector3(dir[1], dir[2], dir[3])$normalize()

  plane_nromals <- new_vector3(
    plane$normals[1, ],
    plane$normals[2, ],
    plane$normals[3, ]
  )
  dot <- plane_nromals$dot(vec3)
  normals <- plane$normals[1:3, , drop = FALSE] * ifelse(dot >= 0, 1, -1)

  raycaster <- vcg_raycaster(
    target, plane$vb[1:3, , drop = FALSE], ray_direction = normals, both_sides = FALSE)
  search_results <- vcg_kdtree_nearest(target, query = plane, k = 1)

  raycaster_valid <- raycaster$has_intersection &
    abs(raycaster$distance) < search_results$distance[, 1] * 2

  raycaster_invalid <- raycaster$has_intersection & !raycaster_valid
  # raycaster$distance[raycaster_invalid]
  # search_results$distance[raycaster_invalid, 1]
  # plane$normals[,1]

  if(any(raycaster_invalid)) {

    # ray-casting not working
    # change ray directions so the direction points towards
    # mean of normals and shortest point
    idx <- search_results$index[raycaster_invalid]
    plane_nromals$from_array(
      target$vb[1:3, idx, drop = FALSE] - plane$vb[1:3, raycaster_invalid, drop = FALSE])
    plane_nromals$normalize()
    plane_nromals$from_array(plane_nromals[] + normals[, raycaster_invalid])$normalize()

    raycaster2 <- vcg_raycaster(
      target, plane$vb[1:3, raycaster_invalid, drop = FALSE],
      ray_direction = plane_nromals[], both_sides = FALSE)

    raycaster2_valid <- raycaster2$has_intersection &
      abs(raycaster2$distance) < search_results$distance[raycaster_invalid, 1] * 2

    if(any(raycaster2_valid)) {
      idx <- which(raycaster_invalid)[raycaster2_valid]
      raycaster$intersection[, idx] <- raycaster2$intersection[, raycaster2_valid]
      raycaster_valid[idx] <- TRUE
    }
  }

  # plane_proj <- plane
  plane$vb[, raycaster_valid] <- raycaster$intersection[, raycaster_valid]
  plane <- vcg_update_normals(plane, weight = "area")
  plane <- vcg_smooth_explicit(plane, type = "taubin", iteration = 2)

  if(diagnostic && system.file(package = "ieegio") != '') {
    # for RAVE users only
    ieegio <- asNamespace('ieegio')
    merged <- merge(ieegio$as_ieegio_surface(plane),
                    ieegio$as_ieegio_surface(target),
                    merge_type = "geometry", verbose = FALSE)

    pal <- grDevices::colorRampPalette(
      c("#053061", "#2166ac", "#4393c3", "#92c5de", "#b4b4b4",
        "#f4a582", "#d6604d", "#b2182b", "#67001f"))(255)
    print(plot(
      merge(
        merged,
        ieegio$as_ieegio_surface(
          measurements = data.frame(
            distance_to_surface = c(
              search_results$distance[, 1], rep(0, ncol(target$vb))
            ),
            has_intersection = c(
              raycaster$has_intersection * 3+1, rep(0, ncol(target$vb))
            )
          )
        ),
        verbose = FALSE
      ),
      col = c("#AAAAAAAA", paste0(pal, "AA")),
      name = c("measurements", "distance_to_surface")
    ))
  }

  t(plane$vb[1:3, , drop = FALSE])

}


