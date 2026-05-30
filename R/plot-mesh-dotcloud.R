# Internal: apply user-supplied world-space clipping planes to a 3 x N matrix
# of points. Returns a logical keep mask of length N. `eye_world` is used
# only for planes whose `side` element is 0 (auto camera-facing). A NULL or
# zero-length `clipping_plane` returns an all-TRUE mask. A bare length-5
# numeric vector is accepted as shorthand for a one-plane list.
apply_clipping_planes <- function(points_world, eye_world, clipping_plane) {
  n_pts <- ncol(points_world)
  keep <- rep(TRUE, n_pts)
  if (!length(clipping_plane)) return(keep)
  if (is.numeric(clipping_plane) && !is.list(clipping_plane)) {
    clipping_plane <- list(clipping_plane)
  } else if (!is.list(clipping_plane)) {
    stop("`clipping_plane` must be NULL, a length-5 numeric vector, or a list of such vectors")
  }
  for (cp in clipping_plane) {
    cp <- as.numeric(cp)
    if (length(cp) != 5L || anyNA(cp)) {
      stop("`clipping_plane`: each plane must be a numeric vector of length 5 with no NA")
    }
    nrm <- cp[1:3]
    nlen <- sqrt(sum(nrm * nrm))
    if (!is.finite(nlen) || nlen == 0) {
      stop("`clipping_plane`: the first three elements (normal vector) must not be zero")
    }
    nrm <- nrm / nlen
    d <- cp[[4L]] / nlen
    sd <- cp[[5L]]
    if (!isTRUE(sd %in% c(-1, 0, 1))) {
      stop("`clipping_plane`: the fifth element must be -1, 0, or 1")
    }
    if (sd == 0) {
      signed_eye <- sum(nrm * eye_world) - d
      sd <- if (signed_eye >= 0) 1 else -1
    }
    signed_pts <- as.numeric(nrm %*% points_world) - d
    if (sd > 0) {
      keep <- keep & (signed_pts >= 0)
    } else {
      keep <- keep & (signed_pts <= 0)
    }
  }
  keep
}

# Internal: resolve auto-computed xlim/ylim with a `zoom` magnification that
# is aware of the current device's plot region (`par("pin")`) and the
# requested aspect ratio `asp`. When both `asp` and `pin` are usable, the
# visible data window after `asp`-induced expansion is shrunk by `zoom`, so
# the limiting axis is zoomed by exactly `zoom` and the other axis is
# zoomed less (filling extra margin) while preserving `asp`. When `asp` or
# `pin` is unavailable, falls back to a symmetric shrink by `zoom` around
# the data midpoint. User-supplied (non-NULL) `xlim`/`ylim` are passed
# through unchanged.
zoom_limits <- function(xlim, ylim, x_all, y_all, zoom, asp) {
  x_user <- length(xlim) > 0
  y_user <- length(ylim) > 0
  if (!x_user) xlim <- range(x_all, na.rm = TRUE)
  if (!y_user) ylim <- range(y_all, na.rm = TRUE)
  if (x_user && y_user) return(list(xlim = xlim, ylim = ylim))
  if (!(is.finite(zoom) && zoom > 0) || zoom == 1) {
    return(list(xlim = xlim, ylim = ylim))
  }

  dx <- xlim[2L] - xlim[1L]
  dy <- ylim[2L] - ylim[1L]
  mx <- mean(xlim)
  my <- mean(ylim)

  pin <- tryCatch(graphics::par("pin"), error = function(e) NULL)
  asp_ok <- length(asp) == 1L && is.finite(asp) && asp > 0
  pin_ok <- !is.null(pin) && length(pin) == 2L &&
            all(is.finite(pin)) && all(pin > 0)
  data_ok <- is.finite(dx) && is.finite(dy) && dx > 0 && dy > 0

  if (asp_ok && pin_ok && data_ok) {
    pw <- pin[[1L]]
    ph <- pin[[2L]]
    # Physical inches per x-data-unit when both axes fit (asp y/x ratio):
    #   scale_x = min(pw/dx, ph/(asp*dy))
    # Visible ranges: vis_dx = pw/scale_x, vis_dy = ph/(asp*scale_x).
    # Zooming in by `zoom` scales `scale_x` by `zoom`.
    scale_x <- min(pw / dx, ph / (asp * dy))
    new_dx <- pw / (scale_x * zoom)
    new_dy <- ph / (asp * scale_x * zoom)
  } else {
    new_dx <- dx / zoom
    new_dy <- dy / zoom
  }

  if (!x_user) xlim <- c(mx - new_dx / 2, mx + new_dx / 2)
  if (!y_user) ylim <- c(my - new_dy / 2, my + new_dy / 2)
  list(xlim = xlim, ylim = ylim)
}

#' @title Render one or more meshes as an orthographic dot cloud in base R
#' @description
#' Projects mesh vertices onto a 2D plane using an orthographic camera defined
#' by an eye position, a look-at point, and an up direction, then draws the
#' projected dots with base-R \code{\link[graphics]{plot}}.  Each dot is
#' rendered opaque, but its \code{cex} is modulated by a rim-light weight
#' \eqn{1 - |n \cdot z_{cam}|}: front- and back-facing vertices shrink toward
#' zero size, while grazing-edge vertices keep their full size.  This
#' size-modulated trick gives a rim-light look without paying R's per-point
#' transparency-blending cost.
#'
#' Multiple meshes share a single depth space: all vertices are projected
#' together and sorted globally by depth before a single \code{plot()} call,
#' so the painter's algorithm works correctly across meshes.  Meshes without
#' faces (point clouds) are rendered at full size and are not affected
#' by the \code{side} filter.
#'
#' @param mesh a \code{'mesh3d'} object, or a list of \code{'mesh3d'} objects.
#'   Meshes without a face matrix (\code{$it}) are treated as point clouds.
#'   If \code{mesh$normals} is absent and the mesh has faces, normals are
#'   recomputed with \code{\link{vcg_update_normals}}.
#' @param eye numeric vector of length 3 - camera position in world space.
#' @param lookat numeric vector of length 3 - the world-space point the camera
#'   is looking at.
#' @param up numeric vector of length 3 - a world-space vector indicating
#'   which direction is "up" for the camera; defaults to \code{c(0, 1, 0)}.
#' @param col base color(s) for the dots.  Accepted forms:
#'   \itemize{
#'     \item A single color string - applied to all meshes.
#'     \item A character vector of length 2 or more (not matching any vertex
#'       count) - used as a depth-ordered color gradient applied to all meshes.
#'     \item A character vector of length equal to the number of vertices in a
#'       single mesh - per-vertex colors for that mesh.
#'     \item A list of length equal to the number of meshes - each element is
#'       one of the above, applied to the corresponding mesh.
#'   }
#'   Default \code{"gray30"}.
#' @param pch point character; a scalar or vector/list with one value per mesh,
#'   recycled as necessary.  Default \code{16L}.
#' @param cex point expansion factor; a scalar or vector/list with one value
#'   per mesh, recycled as necessary.  Default \code{0.1}.
#' @param add logical; if \code{TRUE} the dots are added to an existing plot
#'   instead of opening a new one.  Default \code{FALSE}.
#' @param axes logical; whether to draw axes on a new plot.  Ignored when
#'   \code{add = TRUE}.  Default \code{FALSE}.
#' @param asp aspect ratio of the new plot; default \code{1} (equal scaling).
#'   Ignored when \code{add = TRUE}.
#' @param xlim,ylim axis limits for the new plot; \code{NULL} (default) lets R
#'   choose automatically from all meshes' projected vertices.  Ignored when
#'   \code{add = TRUE}.
#' @param zoom positive numeric magnification applied to the auto-computed
#'   axis limits when \code{xlim} or \code{ylim} is \code{NULL}.  Values
#'   greater than \code{1} zoom in, values in \code{(0, 1)} zoom out.
#'   When \code{asp} is set, the zoom is plot-region aware (queried via
#'   \code{par("pin")}): the axis that already fills the device after
#'   \code{asp}-induced expansion is zoomed by exactly \code{zoom}, and
#'   the other axis is zoomed less so the data fills more of the plot
#'   region while preserving the requested \code{asp}.  Ignored for any
#'   axis whose limit was supplied explicitly, and when \code{add = TRUE}.
#'   Default \code{1}.
#' @param xlab,ylab axis labels for the new plot.  Ignored when
#'   \code{add = TRUE}.
#' @param normal_weight passed to \code{\link{vcg_update_normals}} when normals
#'   must be recomputed; one of \code{"auto"} (area if no normals present,
#'   otherwise skip), \code{"area"}, or \code{"angle"}.  Default \code{"auto"}.
#' @param side which side of meshed surfaces to render.  One of \code{"front"}
#'   (default), \code{"back"}, or \code{"both"} (renders all vertices).  Point clouds are always
#'   rendered regardless of this setting.
#' @param mesh_clipping numeric in \code{[0, 1]} controlling how much of
#'   the surface is kept relative to the camera-facing direction.
#'   Vertices whose rim-light weight (\eqn{1 - |n \cdot z_{cam}|}) is at
#'   or below \code{1 - mesh_clipping} are dropped, leaving only the
#'   silhouette/grazing band.  \code{mesh_clipping = 1} (no clipping)
#'   keeps all vertices; smaller values peel away more of the
#'   front/back-facing dots.  Surviving vertices are drawn opaque, and
#'   their \code{cex} is multiplied by the rim-light weight so
#'   grazing-edge dots appear larger than near-front-facing dots.
#'   Drawing opaque points with size-modulated weight is dramatically
#'   faster than R's true per-point transparency.  Default \code{0.7}.
#' @param alpha numeric in \code{[0, 1]}; one value per mesh (recycled),
#'   sets the transparency of each mesh (\code{1} = fully opaque,
#'   \code{0} = fully transparent).  Vertices belonging to a mesh with
#'   \code{alpha = 0} are dropped entirely.  Default \code{1}.
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
#'   per vertex.  A single length-5 numeric vector is also accepted as
#'   shorthand for a one-plane list.  Default \code{NULL} (no clipping).
#' @param clipping_plane_enabled logical vector, one entry per mesh
#'   (recycled), controlling whether \code{clipping_plane} is applied to
#'   that mesh.  \code{TRUE} (default) means the mesh participates in
#'   clipping; \code{FALSE} exempts the mesh entirely (all of its vertices
#'   are kept regardless of the clipping planes).  Has no effect when
#'   \code{clipping_plane} is \code{NULL}.
#' @param ... additional graphical parameters forwarded to
#'   \code{\link[graphics]{plot.default}} (new plot) or
#'   \code{\link[graphics]{points}} (when \code{add = TRUE}).
#'
#' @return Invisibly returns a list with components \code{xlim} and
#'   \code{ylim} (the plot limits used).
#'
#' @examples
#'
#' mesh <- vcg_isosurface(left_hippocampus_mask)
#'
#' # Side view - rim-light shows the outline of the hippocampus
#' plot_mesh_dotcloud(
#'   mesh,
#'   eye = c(150, 0, 0),
#'   lookat = c(0, 0, 0),
#'   up = c(0, 0, 1),
#'   col = "steelblue",
#'   cex = 2
#' )
#'
#' # Two meshes: surface + electrode point cloud
#' n_elec <- 20
#' electrodes <- structure(
#'   list(vb = matrix(rnorm(3 * n_elec, sd = 5), 3, n_elec) +
#'          rowMeans(mesh$vb)[1:3]),
#'   class = "mesh3d"
#' )
#' plot_mesh_dotcloud(
#'   mesh  = list(mesh, electrodes),
#'   eye    = c(150, 0, 0),
#'   lookat = c(0, 0, 0),
#'   up     = c(0, 0, 1),
#'   col    = list("steelblue", "red"),
#'   pch    = c(16L, 17L),
#'   cex    = c(2, 1.2)
#' )
#'
#' @seealso \code{\link{vcg_update_normals}}, \code{\link{vcg_isosurface}}
#'
#' @inheritSection ensure_mesh3d Coercing Surface Inputs
#'
#' @export
plot_mesh_dotcloud <- function(
    mesh,
    eye = c(0, 0, 1000),
    lookat = c(0, 0, 0),
    up = c(0, 1, 0),
    col = c("white", "gray30"),
    pch = 16L,
    cex = 0.1,
    add = FALSE,
    axes = FALSE,
    asp = 1,
    xlim = NULL,
    ylim = NULL,
    zoom = 1,
    xlab = "",
    ylab = "",
    normal_weight = c("auto", "area", "angle"),
    side = c("front", "back", "both"),
    mesh_clipping = 0.7,
    alpha = 1,
    clipping_plane = NULL,
    clipping_plane_enabled = TRUE,
    ...
) {
  normal_weight <- match.arg(normal_weight)
  side <- match.arg(side)

  # DIPSAUS DEBUG START
  # mesh <- ieegio::read_surface("~/rave_data/raw_dir/yael_demo_001/rave-imaging/fs/surf/lh.pial")
  # normal_weight = "auto"
  # side <- "front"
  # list2env(envir = .GlobalEnv, list(eye = c(0, 0, 1),
  #                                   lookat = c(0, 0, 0),
  #                                   up = c(0, 1, 0),
  #                                   col = "gray30",
  #                                   pch = 16L,
  #                                   cex = 0.1,
  #                                   add = FALSE,
  #                                   axes = FALSE,
  #                                   asp = 1,
  #                                   xlim = NULL,
  #                                   ylim = NULL,
  #                                   xlab = "",
  #                                   ylab = ""))
  # DIPSAUS DEBUG END

  # ---- 0. Normalize mesh to a list --------------------------------------
  if (inherits(mesh, "mesh3d") || (!is.null(mesh$vb))) {
    mesh_list <- list(mesh)
  } else {
    mesh_list <- mesh
  }
  n_meshes <- length(mesh_list)

  # col: list -> per mesh; scalar/vector -> same col spec for every mesh
  if (is.list(col)) {
    col_list <- rep_len(col, n_meshes)
  } else {
    col_list <- rep(list(col), n_meshes)
  }
  # Normalize every color spec to #RRGGBB up front so downstream code never
  # needs col2rgb / adjustcolor again.
  col_list <- lapply(col_list, function(ci) {
    substr(grDevices::adjustcolor(ci), 1L, 7L)
  })

  # pch/cex/alpha: one value per mesh (recycled)
  pch_vec   <- rep_len(as.integer(unlist(pch)),   n_meshes)
  cex_vec   <- rep_len(as.numeric(unlist(cex)),   n_meshes)
  alpha_vec <- pmax(0, pmin(1, rep_len(as.numeric(unlist(alpha)), n_meshes)))
  clip_enabled_vec <- rep_len(as.logical(clipping_plane_enabled), n_meshes)
  clip_enabled_vec[is.na(clip_enabled_vec)] <- TRUE

  # ---- 1. View matrix: world -> camera (R6 Vector3 / Matrix4) -----------
  # look_at builds R_c2w (camera-to-world, rotation only, three.js convention).
  # invert() == transpose() for a pure rotation, giving R_w2c.
  # set_position() adds the translation column: t = -R_w2c * eye.
  eye_v3    <- as_vector3(eye[1:3])
  lookat_v3 <- as_vector3(lookat[1:3])
  up_v3     <- as_vector3(up[1:3])

  rot_w2c <- Matrix4$new()$look_at(eye_v3, lookat_v3, up_v3)$invert()
  t_vec   <- Vector3$new()$from_array(eye[1:3])$apply_matrix4(rot_w2c)$multiply_scalar(-1)
  view    <- rot_w2c$clone2()$set_position(t_vec)

  # ---- 2. Process each mesh: extract geometry and normals ---------------
  # Only the per-mesh normal recomputation needs a loop; everything else
  # below is vectorized over the concatenated arrays.
  mesh_data <- vector("list", n_meshes)

  for (i in seq_len(n_meshes)) {
    m   <- ensure_mesh3d(mesh_list[[i]])
    vb  <- m$vb
    if (nrow(vb) > 3L) vb <- vb[1:3, , drop = FALSE]
    n_v <- ncol(vb)

    has_faces <- !is.null(m$it) && length(m$it) > 0L

    if (has_faces) {
      nw <- if (normal_weight == "auto") "area" else normal_weight
      if (is.null(m$normals) || normal_weight != "auto") {
        m <- vcg_update_normals(m, weight = nw)
      }
      nm <- m$normals
      if (nrow(nm) > 3L) nm <- nm[1:3, , drop = FALSE]
    } else {
      # Placeholder normals for point clouds (zero columns - facing and rim
      # weight are overridden via the has_faces_per_vert mask below).
      nm <- matrix(0, nrow = 3L, ncol = n_v)
    }

    mesh_data[[i]] <- list(
      vb        = vb,
      normals   = nm,
      has_faces = has_faces,
      n_verts   = n_v,
      pch       = pch_vec[[i]],
      cex       = cex_vec[[i]],
      col       = col_list[[i]]
    )
  }

  # ---- 3. Concatenate everything, then run vectorized ops ---------------
  all_vb      <- do.call(cbind, lapply(mesh_data, `[[`, "vb"))      # 3 x N
  all_normals <- do.call(cbind, lapply(mesh_data, `[[`, "normals")) # 3 x N
  n_verts_vec <- vapply(mesh_data, `[[`, 0L, "n_verts")
  total_n     <- sum(n_verts_vec)

  # per-vertex has_faces mask (TRUE for vertices belonging to a meshed surface)
  has_faces_per_vert <- rep(
    vapply(mesh_data, `[[`, logical(1), "has_faces"),
    n_verts_vec
  )

  # Single view transform for ALL vertices
  cam_pts <- Vector3$new()$from_array(all_vb)$apply_matrix4(view)
  x_all     <- cam_pts[1L, ]
  y_all     <- cam_pts[2L, ]
  depth_all <- -cam_pts[3L, ]    # positive = in front of camera

  # Single direction transform for ALL normals (zero columns for point
  # clouds are fine - they get overridden by the has_faces mask).
  cam_normals <- Vector3$new()$from_array(all_normals)$transform_direction(view)
  facing_all  <- cam_normals[3L, ]   # cam Z = -fwd component

  # Rim-light weight: 1 at grazing edges, 0 at front/back.
  # Point clouds: rim_weight = 1 (override) so they render at full size.
  rim_weight <- pmax(0, pmin(1, 1 - abs(facing_all)))
  rim_weight[!has_faces_per_vert] <- 1.0
  facing_all[!has_faces_per_vert] <- NA_real_

  # pch / cex / alpha: replicate the per-mesh scalar n_v times
  pch_all   <- rep(pch_vec,   n_verts_vec)
  cex_all   <- rep(cex_vec,   n_verts_vec)
  alpha_all <- rep(alpha_vec, n_verts_vec)

  # ---- 4. Per-mesh base colors (only the color step needs per-mesh logic
  # because length(col_i) is compared to that mesh's n_v) -----------------
  # Normalized depth used for gradient colors (global across all meshes).
  dr      <- range(depth_all, na.rm = TRUE)
  dr_span <- dr[2L] - dr[1L]
  if (dr_span < .Machine$double.eps) {
    depth_norm_all <- rep(0.5, total_n)
  } else {
    depth_norm_all <- (depth_all - dr[1L]) / dr_span
  }

  col_all <- character(total_n)
  offsets <- cumsum(c(0L, n_verts_vec))

  for (i in seq_len(n_meshes)) {
    n_v <- n_verts_vec[[i]]
    if (n_v == 0L) next
    idx   <- seq.int(offsets[[i]] + 1L, offsets[[i + 1L]])
    col_i <- col_list[[i]]

    if (length(col_i) == n_v) {
      # Per-vertex colors (already #RRGGBB from upfront normalization)
      base_col <- col_i
    } else if (length(col_i) == 1L) {
      base_col <- rep(col_i[[1L]], n_v)
    } else {
      # Multi-color gradient mapped to global depth of this mesh's vertices.
      # colorRampPalette returns #RRGGBB strings already.
      ramp     <- grDevices::colorRampPalette(rev(col_i))(256L)
      base_col <- ramp[ceiling(depth_norm_all[idx] * 255) + 1L]
    }
    col_all[idx] <- base_col
  }

  # Append per-vertex alpha to the #RRGGBB strings -> #RRGGBBAA.
  # Skip the paste entirely when every mesh is fully opaque so the
  # graphics device can take the no-alpha fast path.
  if (any(alpha_vec < 1)) {
    alpha_hex <- sprintf("%02X", as.integer(round(alpha_all * 255)))
    col_all   <- paste0(col_all, alpha_hex)
  }

  # ---- 5. mesh_clipping cull + side filter (point clouds always pass) ---
  keep <- !is.na(rim_weight) & (rim_weight >= (1 - mesh_clipping)) & (alpha_all > 0)
  if (side == "front") {
    keep <- keep & (is.na(facing_all) | facing_all > 0)
  } else if (side == "back") {
    keep <- keep & (is.na(facing_all) | facing_all < 0)
  }

  # User-supplied world-space clipping planes (per-vertex cull). Meshes
  # with `clipping_plane_enabled = FALSE` are exempted from this cull.
  if (length(clipping_plane)) {
    cp_keep <- apply_clipping_planes(all_vb, eye[1:3], clipping_plane)
    clip_enabled_all <- rep(clip_enabled_vec, n_verts_vec)
    cp_keep[!clip_enabled_all] <- TRUE
    keep <- keep & cp_keep
  }

  # ---- 6. Open canvas + screen limits -----------------------------------
  # Call plot.new() before zoom_limits so par("pin") reflects the true
  # plot-region dimensions of this device (margins applied). Then
  # plot.window() locks in the zoomed limits + asp before drawing.
  if (!add) {
    graphics::plot.new()
  }
  lim <- zoom_limits(xlim, ylim, x_all, y_all, zoom, asp)
  xlim <- lim$xlim
  ylim <- lim$ylim
  if (!add) {
    graphics::plot.window(xlim = xlim, ylim = ylim, asp = asp, ...)
    if (axes) {
      graphics::axis(1)
      graphics::axis(2)
      graphics::box()
    }
    if (nzchar(xlab) || nzchar(ylab)) {
      graphics::title(xlab = xlab, ylab = ylab)
    }
  }

  # ---- 7. Filter, painter's sort (far first), fold rim weight into cex --
  # Dots are drawn at full color; cex is scaled by the per-vertex rim
  # weight so grazing-edge dots remain visually heavier than near
  # front-facing dots without paying R's slow transparency-blending cost.
  x_k   <- x_all[keep]
  y_k   <- y_all[keep]
  col_k <- col_all[keep]
  pch_k <- pch_all[keep]
  cex_k <- cex_all[keep] * rim_weight[keep]

  ord   <- order(depth_all[keep], decreasing = TRUE)
  x_k   <- x_k[ord]
  y_k   <- y_k[ord]
  col_k <- col_k[ord]
  pch_k <- pch_k[ord]
  cex_k <- cex_k[ord]

  # ---- 8. Draw points ---------------------------------------------------
  graphics::points(x_k, y_k, col = col_k, pch = pch_k, cex = cex_k, ...)

  invisible(list(xlim = xlim, ylim = ylim))
}

