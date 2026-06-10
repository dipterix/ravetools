# Internal helper: build a sub-mesh from a set of face indices of orig_mesh.
# Vertex indices are remapped; orig_vertex stores the original column indices.
build_patch_mesh <- function(orig_mesh, face_indices) {
  vb <- orig_mesh$vb           # 4×nv (homogeneous)
  it <- orig_mesh$it           # 3×nf, 1-based
  patch_faces <- it[, face_indices, drop = FALSE]
  orig_v_idx  <- sort(unique(as.integer(patch_faces)))

  remap <- integer(ncol(vb))
  remap[orig_v_idx] <- seq_along(orig_v_idx)
  new_it <- matrix(remap[as.integer(patch_faces)], nrow = 3L)
  new_vb <- vb[, orig_v_idx, drop = FALSE]

  m <- list(vb = new_vb, it = new_it)
  if (!is.null(orig_mesh$normals))
    m$normals <- orig_mesh$normals[, orig_v_idx, drop = FALSE]
  m$orig_vertex <- orig_v_idx
  class(m) <- c("ravetools_mesh3d", "mesh3d")
  m
}

#' @title Split a mesh into two patches along a geodesic boundary
#' @description
#' Connects a set of surface \verb{waypoints} with geodesic paths to form a
#' closed boundary loop, then splits the mesh into the two regions created by
#' that loop.
#'
#' Consecutive \verb{waypoints} are joined by the geodesic (Dijkstra) shortest
#' path along the mesh surface; the last \verb{waypoint} connects back to the
#' first.
#'
#' @param mesh triangular mesh of class \code{'mesh3d'}.
#' @param waypoints numeric matrix with exactly 3 columns (\code{x}, \code{y},
#'   \code{z}) and at least 3 rows; each row is a 3-D coordinate defining a
#'   boundary corner. After any mesh subdivision, each coordinate is snapped to
#'   the nearest vertex in the (possibly refined) mesh. Consecutive snapped
#'   vertices are joined by geodesic paths to form the closed boundary loop.
#' @param seed_vertex integer (optional, 1-based). A vertex known to be inside
#'   the desired first patch. When \code{NULL} (default) the smaller of the
#'   two regions is returned first.
#' @param max_edge_length numeric (optional). When positive and finite, the mesh
#'   is refined before patching so that no edge exceeds this length. Global
#'   edge subdivision (\code{\link{vcg_subdivision}}) is applied repeatedly
#'   until the average edge length falls below the threshold (fast), then
#'   \code{\link{vcg_subdivide_max_edge_length}} handles any remaining outlier
#'   edges. \verb{Waypoint} coordinates are snapped to vertices after
#'   refinement, so finer vertices improve boundary accuracy.
#'   Default \code{NA} (no refinement).
#'
#' @returns A length-2 list of \code{mesh3d} objects. Each contains:
#' \describe{
#'   \item{\code{$orig_vertex}}{1-based integer vector: new vertex index \code{i}
#'     corresponds to column \code{orig_vertex[i]} of the original \code{mesh$vb}.}
#' }
#' The first element is the patch whose centroid is closest to the mean
#' \verb{waypoint} position; the second is the complementary connected patch.
#' On a multi-manifold mesh, disconnected components not adjacent to the
#' boundary loop appear in neither patch. When the boundary loop does not
#' divide the mesh (degenerate \verb{waypoints}), the second element is
#' \code{NULL}.
#'
#' @note All \verb{waypoints} must lie on the same connected component of the mesh.
#'   The mesh should be manifold; run \code{\link{vcg_fix_defects}} first
#'   if needed.
#'
#' @seealso \code{\link{dijkstras_surface_distance}}, \code{\link{surface_path}},
#'   \code{\link{vcg_fix_defects}}
#'
#' @export
vcg_mesh_patch <- function(mesh, waypoints, seed_vertex = NULL,
                            max_edge_length = NA) {
  mesh <- meshintegrity(mesh, facecheck = TRUE)

  waypoints <- as.matrix(waypoints)
  if (!is.numeric(waypoints) || ncol(waypoints) != 3L || nrow(waypoints) < 3L)
    stop("vcg_mesh_patch: `waypoints` must be a numeric matrix with 3 columns (x,y,z) and at least 3 rows")

  # ---- optional mesh refinement ----
  if (isTRUE(max_edge_length > 0) && is.finite(max_edge_length)) {
    avg_len <- vcg_average_edge_length(mesh)
    while (avg_len >= max_edge_length) {
      mesh    <- vcg_subdivision(mesh, "edge")
      avg_len <- vcg_average_edge_length(mesh)
    }
    mesh <- vcg_subdivide_max_edge_length(mesh, max_edge_len = max_edge_length)
  }

  vb <- mesh$vb[1:3, , drop = FALSE]
  it <- mesh$it          # 3×nf, 1-based

  # ---- snap waypoint coordinates to nearest vertices (after subdivision) ----
  vb3 <- vb  # 3×nv
  waypoint_verts <- apply(waypoints, 1L, function(coord) {
    which.min(colSums((vb3 - coord)^2))
  })
  if (anyDuplicated(waypoint_verts))
    warning("vcg_mesh_patch: two or more waypoints snapped to the same vertex; boundary may be degenerate")

  # dijkstras_surface_distance expects positions as N×3 (rows = vertices)
  # and faces as M×3 (rows = faces), 1-indexed
  positions <- t(vb)
  faces     <- t(it)

  # ---- build boundary loop via consecutive geodesic segments ----
  n_wp <- length(waypoint_verts)
  boundary_vertices <- integer(0)
  for (i in seq_len(n_wp)) {
    src <- waypoint_verts[i]
    dst <- waypoint_verts[if (i == n_wp) 1L else i + 1L]
    dist_res <- dijkstras_surface_distance(
      positions        = positions,
      faces            = faces,
      start_node       = src,
      face_index_start = 1L
    )
    seg <- surface_path(dist_res, target_node = dst)
    if (seg$path[1L] != src)
      stop(sprintf(
        "vcg_mesh_patch: cannot find geodesic path from waypoint %d (vertex %d) to waypoint %d (vertex %d). Are they on the same connected component?",
        i, src,
        if (i == n_wp) 1L else i + 1L, dst
      ))
    # Drop the last vertex (= first vertex of the next segment) to avoid duplicates
    boundary_vertices <- c(boundary_vertices, seg$path[-length(seg$path)])
  }

  # ---- resolve seed face (0-based for C++) ----
  seed_face <- -1L
  if (!is.null(seed_vertex)) {
    sv <- as.integer(seed_vertex)
    fc <- which(colSums(it == sv) > 0L)
    if (!length(fc))
      stop("vcg_mesh_patch: seed_vertex ", sv, " does not belong to any face")
    seed_face <- fc[1L] - 1L   # convert to 0-based
  }

  # ---- C++ BFS flood fill for fi1 ----
  fi1 <- vcgMeshPatchFaces(
    vb_          = vb,
    it_          = it - 1L,                  # 0-based faces
    boundary_seq = boundary_vertices - 1L,   # 0-based vertices
    seed_face    = seed_face
  )

  # ---- fi2: BFS from the other side of the first boundary edge ----
  # Using setdiff would incorrectly include disconnected components on a
  # multi-manifold mesh. Instead, seed the second BFS from a face adjacent
  # to the first boundary edge that is not in fi1.
  v0_1 <- boundary_vertices[1L]   # already 1-based
  v1_1 <- boundary_vertices[2L]
  seed2_candidates <- setdiff(which(colSums(it == v0_1) > 0L & colSums(it == v1_1) > 0L), fi1)
  fi2 <- if (length(seed2_candidates) == 0L) {
    integer(0)   # boundary lies on mesh boundary, or loop fully encloses the mesh
  } else {
    vcgMeshPatchFaces(
      vb_          = vb,
      it_          = it - 1L,
      boundary_seq = boundary_vertices - 1L,
      seed_face    = seed2_candidates[1L] - 1L   # 0-based
    )
  }

  # ---- construct sub-meshes ----
  p1 <- build_patch_mesh(mesh, fi1)
  p2 <- if (length(fi2) > 0L) build_patch_mesh(mesh, fi2) else NULL

  # ---- order: patch closest to mean waypoint coordinate comes first ----
  mean_wp  <- colMeans(waypoints)   # 3-vector from the n×3 input matrix
  centroid <- function(m3d) rowMeans(m3d$vb[1:3, , drop = FALSE])
  d1 <- sqrt(sum((centroid(p1) - mean_wp)^2))
  if (!is.null(p2)) {
    d2 <- sqrt(sum((centroid(p2) - mean_wp)^2))
    if (d2 < d1) { tmp <- p1; p1 <- p2; p2 <- tmp }
  }

  list(p1, p2)
}
