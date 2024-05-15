#' @name dijkstras-path
#' @title Calculate distances along a surface
#' @description
#' Calculate surface distances of graph or mesh using 'Dijkstra' method.
#' @param start_node integer, row index of \code{positions} on where to start
#' calculating the distances. This integer must be 1-indexed and cannot exceed
#' the total number of \code{positions} rows
#' @param positions numeric matrix with no \code{NA} values. The number of row
#' is the total count of nodes (vertices), and the number of columns represent
#' the node dimension. Each row represents a node.
#' @param faces integer matrix with each row containing indices of nodes. For
#' graphs, \code{faces} is a matrix with two columns defining the connecting
#' edges; for '3D' mesh, \code{faces} is a three-column matrix defining the
#' face index of mesh triangles.
#' @param face_index_start integer, the start of the nodes in \code{faces};
#' please specify this input explicitly if the first node is not contained
#' in \code{faces}.
#' Default is \code{NA} (determined by the minimal number in \code{faces}).
#' The reason to set this input is because some programs use \code{1} to
#' represent the first node, some start from \code{0}.
#' @param max_search_distance numeric, maximum distance to iterate;
#' default is \code{NA},
#' that is to iterate and search the whole mesh
#' @param max_edge_length numeric, maximum edge length to consider;
#' default is \code{NA},
#' that is to consider all edges. In triangular mesh objects, it is possible
#' that some edges are too long to consider as "valid" path. User this input
#' to ignore those edges.
#' @param verbose whether to verbose the progress; default is false
#' @param x distance calculation results returned by
#' \code{dijkstras_surface_distance} function
#' @param target_node the target node number to reach (from the starting node);
#' \code{target_node} is always 1-indexed.
#' @returns \code{dijkstras_surface_distance} returns a list distance
#' table with the meta configurations. \code{surface_path} returns a data frame
#' of the node ID (from \code{start_node} to \code{target_node}) and cumulative
#' distance along the shortest path.
#'
#' @examples
#'
#' # ---- Toy example --------------------
#'
#' # Position is 2D, total 6 points
#' positions <- matrix(runif(6 * 2), ncol = 2)
#'
#' # edges defines connected nodes
#' edges <- matrix(ncol = 2, byrow = TRUE, data = c(
#'   1,2,
#'   2,3,
#'   1,3,
#'   2,4,
#'   3,4,
#'   2,5,
#'   4,5,
#'   2,5,
#'   4,6,
#'   5,6
#' ))
#'
#' # calculate distances
#' ret <- dijkstras_surface_distance(
#'   start_node = 1,
#'   positions = positions,
#'   faces = edges,
#'   face_index_start = 1
#' )
#'
#' # get shortest path from the first node to the last
#' path <- surface_path(ret, target_node = 6)
#'
#' # plot the results
#' from_node <- path$path[-nrow(path)]
#' to_node <- path$path[-1]
#' plot(positions, pch = 16, axes = FALSE,
#'      xlab = "X", ylab = "Y", main = "Dijkstra's shortest path")
#' segments(
#'   x0 = positions[edges[,1],1], y0 = positions[edges[,1],2],
#'   x1 = positions[edges[,2],1], y1 = positions[edges[,2],2]
#' )
#'
#' points(positions[path$path,], col = "steelblue", pch = 16)
#' arrows(
#'   x0 = positions[from_node,1], y0 = positions[from_node,2],
#'   x1 = positions[to_node,1], y1 = positions[to_node,2],
#'   col = "steelblue", lwd = 2, length = 0.1, lty = 2
#' )
#'
#' points(positions[1,,drop=FALSE], pch = 16, col = "orangered")
#' points(positions[6,,drop=FALSE], pch = 16, col = "purple3")
#'
#' # ---- Example with mesh ------------------------------------
#'
#' \dontrun{
#'
#'   # Please install the down-stream package `threeBrain`
#'   # and call library(threeBrain)
#'   # the following code set up the files
#'
#'   read.fs.surface <- internal_rave_function(
#'     "read.fs.surface", "threeBrain")
#'   default_template_directory <- internal_rave_function(
#'     "default_template_directory", "threeBrain")
#'   surface_path <- file.path(default_template_directory(),
#'                             "N27", "surf", "lh.pial")
#'   if(!file.exists(surface_path)) {
#'     internal_rave_function(
#'       "download_N27", "threeBrain")()
#'   }
#'
#'   # Example starts from here --->
#'   # Load the mesh
#'   mesh <- read.fs.surface(surface_path)
#'
#'   # Calculate the path with maximum radius 100
#'   ret <- dijkstras_surface_distance(
#'     start_node = 1,
#'     positions = mesh$vertices,
#'     faces = mesh$faces,
#'     max_search_distance = 100,
#'     verbose = TRUE
#'   )
#'
#'   # get shortest path from the first node to node 43144
#'   path <- surface_path(ret, target_node = 43144)
#'
#'   # plot
#'   from_nodes <- path$path[-nrow(path)]
#'   to_nodes <- path$path[-1]
#'   # calculate colors
#'   pal <- colorRampPalette(
#'     colors = c("red", "orange", "orange3", "purple3", "purple4")
#'   )(1001)
#'   col <- pal[ceiling(
#'     path$distance / max(path$distance, na.rm = TRUE) * 1000
#'   ) + 1]
#'   oldpar <- par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))
#'   for(xdim in c(1, 2, 3)) {
#'     if( xdim < 3 ) {
#'       ydim <- xdim + 1
#'     } else {
#'       ydim <- 3
#'       xdim <- 1
#'     }
#'     plot(
#'       mesh$vertices[, xdim], mesh$vertices[, ydim],
#'       pch = ".", col = "#BEBEBE33", axes = FALSE,
#'       xlab = "P - A", ylab = "S - I", asp = 1
#'     )
#'     segments(
#'       x0 = mesh$vertices[from_nodes, xdim],
#'       y0 = mesh$vertices[from_nodes, ydim],
#'       x1 = mesh$vertices[to_nodes, xdim],
#'       y1 = mesh$vertices[to_nodes, ydim],
#'       col = col
#'     )
#'   }
#'
#'   # plot distance map
#'   distances <- ret$paths$distance
#'   col <- pal[ceiling(distances / max(distances, na.rm = TRUE) * 1000) + 1]
#'   selection <- !is.na(distances)
#'
#'   plot(
#'     mesh$vertices[, 2], mesh$vertices[, 3],
#'     pch = ".", col = "#BEBEBE33", axes = FALSE,
#'     xlab = "P - A", ylab = "S - I", asp = 1
#'   )
#'   points(
#'     mesh$vertices[selection, c(2, 3)],
#'     col = col[selection],
#'     pch = "."
#'   )
#'
#'   # reset graphic state
#'   par(oldpar)
#'
#' }
#'
#'
#'
#'
#'
#'
#' @export
dijkstras_surface_distance <- function(start_node, positions, faces, face_index_start = NA, max_search_distance = NA, max_edge_length = NA, verbose = FALSE) {
  stopifnot(is.matrix(positions))
  stopifnot(is.matrix(faces))

  if(!is.integer(faces)) {
    faces <- array(as.integer(faces), dim = dim(faces))
  }

  n_points <- nrow(positions)
  n_indices <- nrow(faces)
  max_index <- max(faces)
  min_index <- min(faces)

  verbose <- isTRUE(as.logical(verbose))


  if(!isTRUE(max_search_distance > 0 && is.finite(max_search_distance))) {
    max_search_distance <- 0.0
  }
  if(!isTRUE(max_edge_length > 0 && is.finite(max_edge_length))) {
    max_edge_length <- 0.0
  }

  if(anyNA(positions)) {
    stop("Cannot handle NA vertex positions")
  }
  if(is.na(min_index)) {
    stop("Cannot handle NA face index")
  }
  start_node <- as.integer(start_node)
  if(!isTRUE( start_node >= 1 & start_node <= n_points)) {
    stop("`surface_distance_dijkstras`: start_node must be an integer indicating the vertex node from which that distance calculation starts. The acceptable value must be from 1 to the number of nodes")
  }
  # Check the face index start (0 or 1)
  face_index_start <- as.integer(face_index_start)
  if(length(face_index_start) != 1 || is.na(face_index_start)) {
    face_index_start <- min_index
  }

  faces <- faces - face_index_start

  if( max_index - face_index_start + 1 > n_points ) {
    warning("`surface_distance_dijkstras`: face index starts from ", face_index_start,
         ". Found face index [", max_index, "] that is greater than maximum allowed [",
         n_points - 1 + face_index_start, "]. Some mesh triangles will be ignored.")
    max_fidx <- n_points - 1
    sel <- rowSums(faces < 0 | faces > max_fidx) == 0
    faces <- faces[sel, , drop = FALSE]
    n_indices <- nrow(faces)
  }

  positions <- t(positions)
  faces <- t(faces)

  re <- dijkstras_path(
    position = positions,
    index = faces,
    indexOrder = order(faces) - 1L,
    nPoints = n_points,
    nIndices = n_indices,
    startIndex = start_node - 1L,
    maxDistance = max_search_distance,
    maxEdgeLen = max_edge_length,
    verbose = verbose
  )


  # previous_node <- re[[1]] + 1
  distance <- re[[2]]
  distance[distance < 0] <- NA_real_

  structure(
    list(
      paths = data.frame(
        node_id = seq_len(n_points),
        prev_id = re[[1]] + 1,
        distance = distance
      ),
      start_node = start_node,
      max_search_distance = ifelse(max_search_distance <= 0, Inf, max_search_distance),
      max_edge_length = ifelse(max_edge_length <= 0, Inf, max_edge_length),
      n_nodes = n_points,
      n_faces = n_indices
    ),
    class = c("ravetools-surf-dist-dijkstras", "ravetools-surf-dist")
  )
}


# x <- dijkstras_surface_distance(1, mesh$vertices, mesh$faces, max_search_distance = NA, verbose = TRUE)

#' @rdname dijkstras-path
#' @export
surface_path <- function(x, target_node) {

  if(!inherits(x, "ravetools-surf-dist")) {
    stop("`surface_path`: input `x` must be a surface distance result. Please check `?dijkstras_surface_distance` to calculate surface distance first.")
  }

  target_node <- as.integer(target_node)
  if( !isTRUE( target_node >= 1 & target_node <= x$n_nodes )) {
    stop("`dijkstras_findpath`: `target_node` must be 1~", x$n_nodes)
  }

  new_idx <- target_node
  idx <- NULL
  while(!is.na(new_idx)) {
    idx <- c(new_idx, idx)
    new_idx <- x$paths$prev_id[[new_idx]]
  }

  data.frame(
    path = idx,
    distance = x$paths$distance[idx]
  )
}


#' Get external function from 'RAVE'
#' @description
#' Internal function used for examples relative to 'RAVE' project and should
#' not be used directly.
#' @param name function or variable name
#' @param pkg 'RAVE' package name
#' @param inherit passed to \code{\link{get0}}
#' @param on_missing default value to return of no function is found
#' @returns Function object if found, otherwise \code{on_missing}.
#' @export
internal_rave_function <- function(name, pkg, inherit = TRUE, on_missing = NULL) {
  if(!pkg %in% c("raveio", "ravedash", "ravebuiltins", "rave", "threeBrain",
                 "dipsaus", "filearray", "readNSx", "rpymat", "rpyANTs")) {
    stop("`extern_function`: Package [", pkg, "] is not a RAVE package.")
  }
  if(system.file(package = pkg) == "") { return(on_missing) }
  ns <- asNamespace(pkg)
  get0(x = name, envir = ns, inherits = inherit, ifnotfound = on_missing)
}
