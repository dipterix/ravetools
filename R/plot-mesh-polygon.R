#' @title Render one or more meshes as flat-shaded triangles in base R
#' @description
#' Projects each triangular face onto a 2D plane using an orthographic camera,
#' shades it with a single color proportional to how directly it faces the
#' camera (Lambert), depth-sorts all faces across all meshes, and draws them
#' in a single \code{\link[graphics]{polygon}} call.
#'
#' Meshes without faces (point clouds) are substituted by a small sphere
#' (\code{\link{vcg_sphere}}) centered at each point and scaled by
#' \code{cex}; they then participate in the same rendering pipeline as ordinary
#' faced meshes.
#'
#' A camera-facing clipping pass discards triangles whose outward normal
#' points along the camera ray (signed \eqn{n \cdot z_{cam}} \eqn{> 1 - \mathrm{mesh\_clipping}}),
#' peeling the front cap off the surface so the back wall (and any
#' interior meshes) become visible. Set \code{mesh_clipping = 1} to disable
#' clipping. Point-cloud meshes (those rendered as substitute
#' \code{\link{vcg_sphere}} instances) are exempt from this clip so they
#' remain solid even when the enclosing surface is peeled.
#'
#' Multiple meshes share a single depth space: all faces are projected,
#' sorted, and drawn together so the painter's algorithm works correctly
#' across meshes.
#'
#' @param mesh a \code{'mesh3d'} object, or a list of \code{'mesh3d'} objects.
#'   Meshes without a face matrix (\code{$it}) are rendered as sphere
#'   instances (one \code{\link{vcg_sphere}} per vertex, radius \code{cex}).
#' @param eye numeric vector of length 3 - camera position in world space.
#' @param lookat numeric vector of length 3 - the world-space point the camera
#'   is looking at.
#' @param up numeric vector of length 3 - world-space "up" direction; defaults
#'   to \code{c(0, 1, 0)}.
#' @param col base color(s) per mesh.  Same forms as
#'   \code{\link{plot_mesh_dotcloud}}: a single color, a depth-gradient vector,
#'   a per-vertex character vector, or a list of any of these (one element per
#'   mesh).  Default \code{c("white", "gray30")}.
#' @param cex radius of the substitute sphere used for point-cloud meshes
#'   (world units).  Has no effect on meshes that already have faces.  Default
#'   \code{1}.
#' @param sphere_subdivision integer subdivision level forwarded to
#'   \code{\link{vcg_sphere}} when substituting point clouds.  Higher values
#'   give smoother spheres at the cost of more triangles.  Default \code{1}
#'   (80 triangles per sphere).
#' @param add logical; if \code{TRUE} the faces are added to an existing plot
#'   instead of opening a new one.  Default \code{FALSE}.
#' @param axes,asp,xlim,ylim,xlab,ylab passed to \code{\link[graphics]{plot.default}}
#'   when \code{add = FALSE}.
#' @param side which side of each triangle to render.  One of \code{"front"}
#'   (default), \code{"back"}, or \code{"both"} (shows all triangles).
#' @param mesh_clipping numeric in \code{[0, 1]} controlling camera-facing
#'   clipping: any triangle whose signed Lambert dot product
#'   (\eqn{n \cdot z_{cam}}, i.e. positive for front-facing) is at or
#'   above \code{mesh_clipping} is discarded.  A value of \code{1}
#'   (default) keeps the full surface; \code{0.5} peels a 60-degree cap
#'   off the front and reveals the back wall (whose absolute Lambert
#'   shade is still large, so it renders brightly); \code{0} clips
#'   every front-facing triangle and shows the interior only.
#'   Back-facing triangles are never clipped by this rule.  Default
#'   \code{1}.
#' @param alpha numeric in \code{[0, 1]}; one value per mesh (recycled),
#'   sets the transparency of each mesh (\code{1} = fully opaque,
#'   \code{0} = fully transparent).  Faces belonging to a mesh with
#'   \code{alpha = 0} are dropped entirely.  Because faces are drawn in
#'   back-to-front (painter's) order, alpha blending across multiple
#'   meshes follows that order.  Default \code{1}.
#' @param clipping_plane optional list of world-space clipping planes used to
#'   hide parts of the scene.  Each plane is a numeric vector of length 5:
#'   the first three entries are the plane normal \eqn{(n_x, n_y, n_z)}
#'   (must be non-zero; normalized internally), the fourth is the signed
#'   distance from the world origin to the plane along that normal (so the
#'   plane equation is \eqn{n \cdot x = d}), and the fifth indicates which
#'   half-space is kept: \code{1} keeps the front side (the side the normal
#'   points to), \code{-1} keeps the back side, and \code{0} keeps whichever
#'   side currently faces the camera (auto-flipped per call based on
#'   \code{eye}).  Multiple planes are intersected.  Clipping is applied
#'   per face using the face centroid (a face is kept only when its
#'   centroid lies on the kept side of every plane).  A single length-5
#'   numeric vector is also accepted as shorthand for a one-plane list.
#'   Default \code{NULL} (no clipping).
#' @param clipping_plane_enabled logical vector, one entry per mesh
#'   (recycled), controlling whether \code{clipping_plane} is applied to
#'   that mesh.  \code{TRUE} (default) means the mesh participates in
#'   clipping; \code{FALSE} exempts the mesh entirely (all of its faces
#'   are kept regardless of the clipping planes).  Has no effect when
#'   \code{clipping_plane} is \code{NULL}.
#' @param shadow_color color used for fully unlit (grazing/back) faces.
#'   The Lambert shade linearly interpolates from \code{shadow_color}
#'   (at \code{shade = 0}) toward the face color lit by a white light
#'   (at \code{shade = 1}), so unlit faces fade to this color instead
#'   of an implicit black background.  Default \code{par("fg")}, so the
#'   shadow contrasts with the current device's canvas.
#' @param light_intensity non-negative scalar controlling the brightness
#'   of the (white) light source: at \code{shade = 1} the face color is
#'   multiplied by \code{light_intensity} (clamped to \code{[0, 1]}).
#'   \code{1} (default) reproduces the face color exactly when fully
#'   lit; smaller values darken the whole mesh, larger values saturate.
#' @param ambient_intensity scalar in \code{[0, 1]} acting as a lower
#'   bound on the Lambert shade so grazing/back faces still receive some
#'   light and never collapse fully to \code{shadow_color}.  An effective
#'   shade of \code{max(shade, ambient_intensity)} is used in the
#'   background-to-light blend.  Default \code{0.2}.
#' @param ... additional graphical parameters forwarded to
#'   \code{\link[graphics]{plot.default}} (new plot only).
#'
#' @return Invisibly returns a list with components \code{xlim} and
#'   \code{ylim} (the plot limits used).
#'
#' @details
#' Limitations of the base-R polygon path (no \code{rgl}):
#' \itemize{
#'   \item Flat shading only - one color per triangle.  No per-vertex
#'     color interpolation.
#'   \item No depth buffer - faces are depth-sorted by centroid (painter's
#'     algorithm).  Interpenetrating triangles can render in the wrong
#'     order.
#'   \item Anti-aliasing seams can appear between adjacent triangles on
#'     raster devices; Cairo-based devices (\code{png(type = "cairo")},
#'     \code{svg()}, \code{pdf()}) produce cleaner output than the default
#'     quartz/X11 path.
#' }
#'
#' @examples
#'
#' mesh <- vcg_isosurface(left_hippocampus_mask)
#'
#' # Surface alone
#' plot_mesh_polygon(
#'   mesh,
#'   eye    = c(150, 0, 0),
#'   lookat = c(0, 0, 0),
#'   up     = c(0, 0, 1),
#'   col    = "steelblue"
#' )
#'
#' # Surface + electrode point cloud (rendered as small icospheres)
#' n_elec <- 20
#' electrodes <- structure(
#'   list(vb = matrix(rnorm(3 * n_elec, sd = 5), 3, n_elec)),
#'   class = "mesh3d"
#' )
#' plot_mesh_polygon(
#'   mesh = list(mesh, electrodes),
#'   eye    = c(150, 0, 0),
#'   lookat = c(0, 0, 0),
#'   up     = c(0, 0, 1),
#'   col    = list("steelblue", "red"),
#'   cex    = 1.5
#' )
#'
#' @seealso \code{\link{plot_mesh_dotcloud}}, \code{\link{vcg_isosurface}}
#'
#' @inheritSection ensure_mesh3d Coercing \verb{ieegio_surface} inputs
#'
#' @export
plot_mesh_polygon <- function(
    mesh,
    eye    = c(0, 0, 1000),
    lookat = c(0, 0, 0),
    up     = c(0, 1, 0),
    col    = c("white", "gray30"),
    cex    = 1,
    add    = FALSE,
    axes   = FALSE,
    asp    = 1,
    xlim   = NULL,
    ylim   = NULL,
    xlab   = "",
    ylab   = "",
    side   = c("front", "back", "both"),
    mesh_clipping = 1,
    sphere_subdivision = 1L,
    alpha  = 1,
    shadow_color = NULL,
    light_intensity = 1,
    ambient_intensity = 0.2,
    clipping_plane = NULL,
    clipping_plane_enabled = TRUE,
    ...
) {
  side <- match.arg(side)

  # DIPSAUS DEBUG START
  # mesh <- vcg_isosurface(ravetools::left_hippocampus_mask)
  # list2env(envir = .GlobalEnv, list(
  #   eye = c(150, 0, 0), lookat = c(0, 0, 0), up = c(0, 0, 1),
  #   col = "steelblue", cex = 1, add = FALSE, axes = FALSE, asp = 1,
  #   xlim = NULL, ylim = NULL, xlab = "", ylab = "", mesh_clipping = 1
  # ))
  # DIPSAUS DEBUG END

  # ---- 0. Normalize mesh + col list -------------------------------------
  if (inherits(mesh, "mesh3d") || (!is.null(mesh$vb))) {
    mesh_list <- list(mesh)
  } else {
    mesh_list <- mesh
  }
  n_meshes <- length(mesh_list)

  if (is.list(col)) {
    col_list <- rep_len(col, n_meshes)
  } else {
    col_list <- rep(list(col), n_meshes)
  }
  col_list <- lapply(col_list, function(ci) {
    substr(grDevices::adjustcolor(ci), 1L, 7L)
  })

  alpha_vec <- pmax(0, pmin(1, rep_len(as.numeric(unlist(alpha)), n_meshes)))
  clip_enabled_vec <- rep_len(as.logical(clipping_plane_enabled), n_meshes)
  clip_enabled_vec[is.na(clip_enabled_vec)] <- TRUE

  # Shading parameters: shadow color (defaults to par("fg")) and white-light
  # intensity. Resolved once here so par() is only queried per call.
  if (is.null(shadow_color)) {
    shadow_color <- graphics::par("fg")
  }
  shadow_rgb      <- grDevices::col2rgb(shadow_color)[, 1L] / 255
  light_intensity <- max(0, as.numeric(light_intensity)[[1L]])
  ambient_intensity <- max(0, min(1, as.numeric(ambient_intensity)[[1L]]))

  cex <- as.numeric(cex)[[1L]]
  sphere_subdivision <- as.integer(sphere_subdivision)[[1L]]

  # ---- 1. View matrix (R6, same as plot_mesh_dotcloud) ------------------
  rot_w2c <- Matrix4$new()$look_at(
    as_vector3(eye[1:3]),
    as_vector3(lookat[1:3]),
    as_vector3(up[1:3])
  )$invert()
  t_vec <- Vector3$new()$from_array(eye[1:3])$apply_matrix4(rot_w2c)$multiply_scalar(-1)
  view  <- rot_w2c$clone2()$set_position(t_vec)

  # ---- 2. Extract vb/faces (substitute small spheres for point clouds) ----
  sphere_tmpl <- vcg_sphere(sphere_subdivision, normals = FALSE)
  sphere_v <- sphere_tmpl$vb[1:3, , drop = FALSE] * cex   # 3 x n_sv (scaled)
  sphere_f <- sphere_tmpl$it                              # 3 x n_sf
  storage.mode(sphere_f) <- "integer"
  n_sv <- ncol(sphere_v)
  n_sf <- ncol(sphere_f)

  mesh_data <- vector("list", n_meshes)
  for (i in seq_len(n_meshes)) {
    m  <- ensure_mesh3d(mesh_list[[i]])
    vb <- m$vb
    if (nrow(vb) > 3L) vb <- vb[1:3, , drop = FALSE]

    has_faces <- !is.null(m$it) && length(m$it) > 0L

    if (has_faces) {
      it <- m$it
      storage.mode(it) <- "integer"
    } else {
      # Instance a vcg_sphere at each point
      n_pts <- ncol(vb)
      if (n_pts == 0L) {
        vb <- matrix(numeric(0), nrow = 3L, ncol = 0L)
        it <- matrix(integer(0), nrow = 3L, ncol = 0L)
      } else {
        # Vertices: repeat sphere verts and add point position
        rep_sv  <- sphere_v[, rep.int(seq_len(n_sv), n_pts), drop = FALSE]
        rep_pts <- vb[, rep(seq_len(n_pts), each = n_sv), drop = FALSE]
        vb <- rep_sv + rep_pts
        # Faces: repeat sphere_f, offset by n_sv per point
        it_flat <- rep.int(as.integer(sphere_f), n_pts) +
          rep(seq.int(0L, by = n_sv, length.out = n_pts), each = 3L * n_sf)
        it <- matrix(it_flat, nrow = 3L)
      }
    }

    mesh_data[[i]] <- list(
      vb         = vb,
      it         = it,
      col        = col_list[[i]],
      n_verts    = ncol(vb),
      n_faces    = ncol(it),
      pointcloud = !has_faces
    )
  }

  # ---- 3. Concatenate vb and faces (offset face indices) ----------------
  n_verts_vec <- vapply(mesh_data, `[[`, 0L, "n_verts")
  n_faces_vec <- vapply(mesh_data, `[[`, 0L, "n_faces")
  pc_vec      <- vapply(mesh_data, `[[`, NA, "pointcloud")
  # Per-face flag: was this face part of an originally-point-cloud mesh?
  is_pointcloud <- rep.int(pc_vec, n_faces_vec)
  # Per-face alpha (one value per mesh, replicated across that mesh's faces)
  alpha_face <- rep(alpha_vec, n_faces_vec)
  total_verts <- sum(n_verts_vec)
  total_faces <- sum(n_faces_vec)

  all_vb <- do.call(cbind, lapply(mesh_data, `[[`, "vb"))     # 3 x V
  vert_offsets <- c(0L, cumsum(n_verts_vec)[-n_meshes])
  if (length(vert_offsets) < n_meshes) vert_offsets <- 0L
  all_it <- do.call(cbind, lapply(seq_len(n_meshes), function(i) {
    if (n_faces_vec[[i]] == 0L) return(matrix(integer(0), 3L, 0L))
    mesh_data[[i]]$it + vert_offsets[[i]]
  }))

  if (total_faces == 0L) {
    if (!add) {
      if (!length(xlim)) xlim <- c(-1, 1)
      if (!length(ylim)) ylim <- c(-1, 1)
      graphics::plot.default(NA, NA, xlim = xlim, ylim = ylim,
                             asp = asp, axes = axes, type = "n",
                             xlab = xlab, ylab = ylab, ...)
    }
    return(invisible(list(xlim = xlim, ylim = ylim)))
  }

  # ---- 4. Single view transform on all vertices -------------------------
  cam_pts <- Vector3$new()$from_array(all_vb)$apply_matrix4(view)
  x_all <- cam_pts[1L, ]
  y_all <- cam_pts[2L, ]
  z_all <- cam_pts[3L, ]

  # ---- 5. Per-face: normal (cam space), facing, centroid depth ----------
  i1 <- all_it[1L, ]
  i2 <- all_it[2L, ]
  i3 <- all_it[3L, ]

  e1x <- x_all[i2] - x_all[i1]
  e1y <- y_all[i2] - y_all[i1]
  e1z <- z_all[i2] - z_all[i1]
  e2x <- x_all[i3] - x_all[i1]
  e2y <- y_all[i3] - y_all[i1]
  e2z <- z_all[i3] - z_all[i1]

  # Cross product, normalize z-component for cos(angle to camera axis)
  nx <- e1y * e2z - e1z * e2y
  ny <- e1z * e2x - e1x * e2z
  nz <- e1x * e2y - e1y * e2x
  n_len <- sqrt(nx * nx + ny * ny + nz * nz)
  facing <- nz / n_len    # in [-1, 1]; NaN for degenerate triangles

  depth_face <- -(z_all[i1] + z_all[i2] + z_all[i3]) / 3   # +ve = in front

  # Lambert shade: |facing| works for both sides (back-faces are covered by
  # front-faces after painter's sort anyway, so this matches a closed mesh).
  shade <- abs(facing)
  shade[is.na(shade)] <- 0

  # ---- 6. Camera-facing clip (peel away the cap facing the camera) -----
  #         + side filter
  # Discard triangles whose outward normal points along the camera ray:
  # signed facing = n . z_cam >= mesh_clipping. mesh_clipping = 1 keeps
  # the full surface; smaller values peel a wider cap and reveal the back
  # wall (whose Lambert shade |facing| is still large, so it renders
  # brightly).
  keep <- facing <= mesh_clipping
  if (side == "front") {
    keep <- keep & (facing > 0)
  } else if (side == "back") {
    keep <- keep & (facing < 0)
  }
  keep <- is_pointcloud | keep   # point clouds are exempt from camera-facing clip
  keep <- keep & (alpha_face > 0)

  # User-supplied world-space clipping planes (per-face cull via centroid).
  # Applies to point-cloud sphere instances as well. Meshes with
  # `clipping_plane_enabled = FALSE` are exempted from this cull.
  if (length(clipping_plane)) {
    centroids <- (all_vb[, i1, drop = FALSE] +
                  all_vb[, i2, drop = FALSE] +
                  all_vb[, i3, drop = FALSE]) / 3
    cp_keep <- .apply_clipping_planes(centroids, eye[1:3], clipping_plane)
    clip_enabled_face <- rep(clip_enabled_vec, n_faces_vec)
    cp_keep[!clip_enabled_face] <- TRUE
    keep <- keep & cp_keep
  }

  # ---- 7. Per-face base color (one entry per mesh) ----------------------
  dr <- range(depth_face, na.rm = TRUE)
  dr_span <- dr[2L] - dr[1L]
  depth_face_norm <- if (dr_span < .Machine$double.eps) {
    rep(0.5, total_faces)
  } else {
    (depth_face - dr[1L]) / dr_span
  }

  col_face       <- character(total_faces)
  face_offsets   <- cumsum(c(0L, n_faces_vec))

  for (i in seq_len(n_meshes)) {
    n_f <- n_faces_vec[[i]]
    if (n_f == 0L) next
    fidx  <- seq.int(face_offsets[[i]] + 1L, face_offsets[[i + 1L]])
    col_i <- col_list[[i]]
    n_v_i <- n_verts_vec[[i]]

    if (length(col_i) == n_v_i) {
      # Per-vertex: take face's first vertex color (flat shading)
      face_v1_local <- all_it[1L, fidx] - vert_offsets[[i]]
      col_face[fidx] <- col_i[face_v1_local]
    } else if (length(col_i) == 1L) {
      col_face[fidx] <- col_i[[1L]]
    } else {
      ramp <- grDevices::colorRampPalette(rev(col_i))(256L)
      col_face[fidx] <- ramp[ceiling(depth_face_norm[fidx] * 255) + 1L]
    }
  }

  # ---- 8. Filter + painter's sort (far first) ---------------------------
  keep_idx <- which(keep)
  if (!length(keep_idx)) {
    if (!add) {
      if (!length(xlim)) xlim <- range(x_all, na.rm = TRUE)
      if (!length(ylim)) ylim <- range(y_all, na.rm = TRUE)
      graphics::plot.default(NA, NA, xlim = xlim, ylim = ylim,
                             asp = asp, axes = axes, type = "n",
                             xlab = xlab, ylab = ylab, ...)
    }
    return(invisible(list(xlim = xlim, ylim = ylim)))
  }

  ord <- keep_idx[order(depth_face[keep_idx], decreasing = TRUE)]

  i1k <- all_it[1L, ord]
  i2k <- all_it[2L, ord]
  i3k <- all_it[3L, ord]
  shade_k <- shade[ord]
  col_k   <- col_face[ord]
  alpha_k <- alpha_face[ord]

  # ---- 9. Apply Lambert shading -----------------------------------------
  # Light is white at brightness `light_intensity`; the fully-lit color
  # is therefore col * light_intensity. The fully-unlit color is
  # `shadow_color` (default par("fg")). The Lambert term `shade_k` is
  # lower-bounded by `ambient_intensity` so rim/back faces still receive
  # some light:
  #   eff_shade = max(shade, ambient_intensity)
  #   shaded    = shadow + eff_shade * (col * light_intensity - shadow)
  # The final pmin clamp is applied to the output channels (not lit_*),
  # so values stay in [0, 1] only after the lerp.
  rgb_mat <- grDevices::col2rgb(col_k) / 255
  lit_r   <- pmin(1, rgb_mat[1L, ] * light_intensity)
  lit_g   <- pmin(1, rgb_mat[2L, ] * light_intensity)
  lit_b   <- pmin(1, rgb_mat[3L, ] * light_intensity)

  eff_shade <- pmax(shade_k, ambient_intensity)

  r_out <- shadow_rgb[1L] + eff_shade * (lit_r - shadow_rgb[1L])
  g_out <- shadow_rgb[2L] + eff_shade * (lit_g - shadow_rgb[2L])
  b_out <- shadow_rgb[3L] + eff_shade * (lit_b - shadow_rgb[3L])

  # Skip the alpha channel when every mesh is fully opaque (faster device path).
  if (any(alpha_vec < 1)) {
    col_shaded <- grDevices::rgb(r_out, g_out, b_out, alpha = alpha_k)
  } else {
    col_shaded <- grDevices::rgb(r_out, g_out, b_out)
  }

  # ---- 10. Plot limits, build NA-separated polygon arrays ---------------
  if (!length(xlim)) xlim <- range(x_all, na.rm = TRUE)
  if (!length(ylim)) ylim <- range(y_all, na.rm = TRUE)

  # 4 rows per triangle: v1, v2, v3, NA. Flatten column-major for polygon().
  tri_x <- rbind(x_all[i1k], x_all[i2k], x_all[i3k], NA_real_)
  tri_y <- rbind(y_all[i1k], y_all[i2k], y_all[i3k], NA_real_)

  # ---- 11. Render -------------------------------------------------------
  if (!add) {
    graphics::plot.default(NA, NA, xlim = xlim, ylim = ylim,
                           asp = asp, axes = axes, type = "n",
                           xlab = xlab, ylab = ylab, ...)
  }
  grDevices::dev.hold()
  on.exit(grDevices::dev.flush(), add = TRUE)
  # Hide the anti-aliasing seam between adjacent triangles by drawing the
  # border in the fill color - but only when fully opaque. With alpha < 1,
  # the border stroke would blend on top of the fill and produce visibly
  # darker triangle edges, so fall back to no border in that case.
  border_col <- if (any(alpha_vec < 1)) NA else col_shaded
  graphics::polygon(c(tri_x), c(tri_y),
                    col = col_shaded, border = border_col)

  invisible(list(xlim = xlim, ylim = ylim))
}
