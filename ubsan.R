# Workload for ubsan.sh -- see that script for how the sanitizer build is made.
#
# Everything runs in ONE R session on purpose. vcglib hands out "user bit"
# flags from function-local statics that live for the whole process, so
# defects of that shape (the CountNonManifoldEdgeFF leak fixed in 0.3.2) only
# appear after several calls and are invisible to any single-call test.
#
# Coverage is AUDITED, not assumed. Section 5 walks the package call graph,
# derives which exported functions can reach a vcglib entry point, and exits
# non-zero if any of them was never exercised -- or was attempted but only
# ever raised an error. Adding a mesh function to ravetools without covering
# it here is therefore a hard failure, not a silent gap. That audit is what
# found plot_mesh_polygon and project_plane, which reach vcglib indirectly
# and would never have been guessed from the vcg_* naming convention.

library(ravetools)

cat("ravetools:", as.character(packageVersion("ravetools")),
    "from", dirname(system.file(package = "ravetools")), "\n")

section <- function(x) cat("\n========== ", x, " ==========\n", sep = "")

pkg <- Sys.getenv("RAVETOOLS_SRC", unset = ".")

attempted <- new.env(parent = emptyenv())
succeeded <- new.env(parent = emptyenv())

# ok() records coverage only when the call actually ran. Silently swallowing an
# error and still claiming coverage would be worse than no audit at all.
ok <- function(name, expr) {
  assign(name, TRUE, envir = attempted)
  good <- TRUE
  tryCatch(suppressWarnings(suppressMessages(force(expr))),
           error = function(e) good <<- FALSE)
  if (good) assign(name, TRUE, envir = succeeded)
  invisible(NULL)
}

# for extra variant calls whose coverage is already established by an ok()
quiet <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(force(expr))),
           error = function(e) invisible(NULL))
}

# ---- 1. full test suite ------------------------------------------------
section("testthat::test_local()")
try(testthat::test_local(path = pkg, stop_on_failure = FALSE,
                         reporter = "summary"))

# ---- 2. the vignette CRAN re-builds ------------------------------------
section("vignettes/mesh-geometry.Rmd")
vig <- file.path(pkg, "vignettes", "mesh-geometry.Rmd")
if (file.exists(vig) && requireNamespace("rmarkdown", quietly = TRUE)) {
  quiet(rmarkdown::render(vig, output_file = file.path(tempdir(), "mg.html"),
                          quiet = TRUE,
                          envir = new.env(parent = globalenv())))
} else {
  cat("skipped (vignette or rmarkdown unavailable)\n")
}

# ---- 3. repeat-call stress loop ----------------------------------------
# Exhaust process-lifetime state. The 0.3.2 bit-flag bug needed four calls to
# corrupt the face allocator, so fewer than ~5 iterations would miss it.
section("repeat-call stress loop")

N     <- 25L   # cheap core operations
N_VAR <- 5L    # expensive / variant-selecting paths (still > 4)

mesh <- vcg_sphere()

# a deliberately defective mesh: duplicating a face makes its three edges
# non-manifold, which drives the repair and diagnostic paths
broken <- mesh
broken$it <- cbind(broken$it, broken$it[, 1, drop = FALSE])

vol <- array(0, c(24, 24, 24))
vol[8:16, 8:16, 8:16] <- 1

# vertices as rows (n x 3), and a logical vertex mask over the whole mesh
verts <- t(mesh$vb[1:3, , drop = FALSE])
vsel <- rep(FALSE, ncol(mesh$vb))
vsel[seq_len(200)] <- TRUE

for (i in seq_len(N)) {
  ok("vcg_sphere",                    vcg_sphere(sub_division = 2))
  ok("vcg_mesh_volume",               vcg_mesh_volume(mesh))
  quiet(vcg_mesh_volume(broken))      # throws; exercises the error path
  ok("vcg_count_edge_defects",        vcg_count_edge_defects(mesh))
  quiet(vcg_count_edge_defects(broken))
  ok("vcg_average_edge_length",       vcg_average_edge_length(mesh))
  ok("vcg_max_edge_length",           vcg_max_edge_length(mesh))
  ok("vcg_fix_defects",               vcg_fix_defects(broken, verbose = FALSE))
  ok("vcg_subdivide_max_edge_length",
     vcg_subdivide_max_edge_length(mesh,
       max_edge_len = vcg_max_edge_length(mesh) * 0.6))
  # selector is a logical mask over ALL vertices, not a vector of indices
  ok("vcg_subset_vertex",             vcg_subset_vertex(mesh, vsel))
  # kdtree takes points as rows (n x 3), not as columns
  ok("vcg_kdtree_nearest",
     vcg_kdtree_nearest(verts, verts[1:10, , drop = FALSE], k = 3))
  ok("vcg_raycaster",
     vcg_raycaster(mesh,
       rbind(mesh$vb[1:3, 1:5, drop = FALSE] * 2, matrix(1, 1, 5)),
       matrix(-1, 3, 5)))
}
cat("core loop: completed", N, "iterations\n")

# ---- 3b. variant paths -------------------------------------------------
# Each argument below selects a DIFFERENT vcglib algorithm. Calling only the
# defaults would leave most of the library unexercised.
section("variant / algorithm-selecting paths")

for (i in seq_len(N_VAR)) {
  # six distinct smoothing algorithms
  for (ty in c("taubin", "laplace", "HClaplace", "fujiLaplace",
               "angWeight", "surfPreserveLaplace")) {
    ok("vcg_smooth_explicit", vcg_smooth_explicit(mesh, type = ty, iteration = 2))
  }
  ok("vcg_smooth_implicit", vcg_smooth_implicit(mesh))
  quiet(vcg_smooth_implicit(mesh, use_mass_matrix = FALSE, fix_border = TRUE,
                            use_cot_weight = TRUE))
  # both subdivision back-ends
  ok("vcg_subdivision", vcg_subdivision(mesh, method = "edge"))
  quiet(vcg_subdivision(mesh, method = "barycenter"))
  # both normal weighting schemes
  ok("vcg_update_normals", vcg_update_normals(mesh, weight = "area"))
  quiet(vcg_update_normals(mesh, weight = "angle"))
  # resampler flag combinations
  ok("vcg_uniform_remesh",
     vcg_uniform_remesh(mesh, voxel_size = 0.3, verbose = FALSE))
  quiet(vcg_uniform_remesh(mesh, voxel_size = 0.3, discretize = TRUE,
                           multi_sample = TRUE, absolute_distance = TRUE,
                           merge_clost = TRUE, verbose = FALSE))
  # isosurface with and without an upper bound
  ok("vcg_isosurface", vcg_isosurface(vol, threshold_lb = 0.5))
  quiet(vcg_isosurface(vol, threshold_lb = 0.5, threshold_ub = 2))
  quiet(vcg_raycaster(mesh,
    rbind(mesh$vb[1:3, 1:5, drop = FALSE] * 2, matrix(1, 1, 5)),
    matrix(-1, 3, 5), both_sides = TRUE))
}
cat("variant loop: completed", N_VAR, "iterations\n")

# ---- 3c. collision detection, all geometry modes -----------------------
# mode_x / mode_y each select a different vcglib traversal (point cloud,
# polyline chain, triangle mesh); mesh-vs-mesh alone covers a third of it.
section("collision detection modes")

pts  <- t(mesh$vb[1:3, seq_len(60), drop = FALSE])
segs <- rbind(pts[1:10, , drop = FALSE], NA, pts[11:20, , drop = FALSE])
geom <- list(list(mesh, "mesh"), list(pts, "points"), list(segs, "segments"))

for (i in seq_len(N_VAR)) {
  for (a in geom) for (b in geom) {
    ok("vcg_detect_collision",
       vcg_detect_collision(a[[1]], b[[1]], mode_x = a[[2]], mode_y = b[[2]]))
    quiet(vcg_detect_collision(a[[1]], b[[1]], mode_x = a[[2]],
                               mode_y = b[[2]], radius = 0.1))
  }
}
cat("collision loop: completed", N_VAR, "iterations x 9 mode pairs\n")

# ---- 3d. geodesic / patch / volume paths -------------------------------
# vcgDijkstra and vcgMeshPatchFaces live in vcgCommon.cpp but are reached from
# dijkstras-path.R and mesh-patch.R, not from vcg.R.
section("geodesic, patch and volume paths")

positions <- matrix(runif(6 * 2), ncol = 2)
edges <- matrix(ncol = 2, byrow = TRUE, data = c(
  1,2, 2,3, 1,3, 2,4, 3,4, 2,5, 4,5, 4,6, 5,6))

for (i in seq_len(N)) {
  ok("dijkstras_surface_distance",
     dijkstras_surface_distance(start_node = 1, positions = positions,
                                faces = edges, face_index_start = 1))
  ret <- quiet(dijkstras_surface_distance(start_node = 1, positions = positions,
                                          faces = edges, face_index_start = 1))
  if (!is.null(ret)) ok("surface_path", surface_path(ret, target_node = 6))

  # same, on a real mesh
  ret2 <- quiet(dijkstras_surface_distance(
    start_node = 1, positions = t(mesh$vb[1:3, , drop = FALSE]),
    faces = t(mesh$it), face_index_start = 1))
  if (!is.null(ret2)) quiet(surface_path(ret2, target_node = 100))

  ok("vcg_mesh_patch", vcg_mesh_patch(mesh, waypoints = diag(1, 3)))
}

for (i in seq_len(N_VAR)) {
  ok("mesh_from_volume",
     mesh_from_volume(vol, threshold = 0.5, verbose = FALSE))
  quiet(mesh_from_volume(vol, threshold = 0.5, verbose = FALSE, remesh = FALSE))
}
cat("geodesic/patch loop: done\n")

# ---- 3e. mris_* surface pipeline ---------------------------------------
# mris.cpp itself contains no vcglib, but mris_remesh reaches it through
# vcg_average_edge_length -- exactly the sort of indirect route the audit
# exists to catch.
section("mris_* surface pipeline")
for (i in seq_len(N_VAR)) {
  ok("mris_curvature", mris_curvature(mesh))
  ok("mris_smooth",    mris_smooth(mesh, n_iter = 2))
  ok("mris_remesh",    mris_remesh(mesh))
  ok("mris_inflate",   mris_inflate(mesh, n_iter = 2))
  ok("mris_sphere",    mris_sphere(mris_inflate(mesh, n_iter = 2), n_iter = 2))
}
cat("mris loop: done\n")

# ---- 3f. rendering paths -----------------------------------------------
# These reach vcglib through ensure_mesh3d / vcg_update_normals / vcg_sphere.
# project_plane chains four entry points on its own.
section("rendering and projection paths")

grDevices::pdf(file = tempfile(fileext = ".pdf"))
on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

for (i in seq_len(N_VAR)) {
  ok("plot_mesh_polygon",  plot_mesh_polygon(mesh, eye = c(10, 10, 10)))
  ok("plot_mesh_dotcloud", plot_mesh_dotcloud(mesh, eye = c(10, 10, 10)))
  ok("project_plane", project_plane(
       target = mesh, width = 8, height = 8, shape = c(16L, 16L),
       initial_positions = t(mesh$vb[1:3, seq_len(256), drop = FALSE]),
       n_iters = 2))
  if (requireNamespace("rgl", quietly = TRUE)) {
    ok("rgl_plot_normals", rgl_plot_normals(mesh))
  }
}
try(grDevices::dev.off(), silent = TRUE)
cat("rendering loop: done\n")

# ---- 4. coverage audit -------------------------------------------------
# Derive, from the installed package itself, which exported functions can
# reach a vcglib entry point, and fail if any was not exercised. This is what
# keeps the script from going stale as ravetools grows.
section("coverage audit")

# The .Call entry points defined in the only two TUs that include vcglib
# (vcgCommon.cpp and vcgCollision.cpp -- verified with `clang++ -MM`).
vcg_entry <- c(
  "vcgAverageEdgeLength", "vcgCountEdgeDefects", "vcgDetectCollision",
  "vcgDijkstra", "vcgEdgeLengthSubdivision", "vcgEdgeSubdivision",
  "vcgFixDefects", "vcgIsoSurface", "vcgKDTreeSearch", "vcgMaxEdgeLength",
  "vcgMeshPatchFaces", "vcgRaycaster", "vcgSmooth", "vcgSmoothImplicit",
  "vcgSphere", "vcgSubset", "vcgUniformResample", "vcgUpdateNormals",
  "vcgVolume"
)

ns <- asNamespace("ravetools")
fns <- Filter(function(n) is.function(get0(n, envir = ns, inherits = FALSE)),
              ls(ns, all.names = TRUE))

# Direct callees of each package function. all.names() on the parsed body
# returns exact symbol names, so this needs no regex -- important because the
# namespace holds operator methods such as `[.Matrix4` whose names are not
# valid patterns.
is_empty_symbol <- function(d) {
  is.symbol(d) && !nzchar(as.character(d))
}
calls_of <- lapply(fns, function(n) {
  fn <- get0(n, envir = ns, inherits = FALSE)
  b <- body(fn)
  syms <- if (is.null(b)) character(0) else all.names(b)
  # formals() yields the empty symbol for arguments with no default; it is a
  # symbol but all.names() cannot be applied to it, so drop it first.
  defaults <- as.list(formals(fn))
  defaults <- defaults[!vapply(defaults, is_empty_symbol, logical(1))]
  for (d in defaults) syms <- c(syms, all.names(d))
  unique(syms)
})
names(calls_of) <- fns

# transitive closure: which functions can reach any vcglib entry point
reaches <- vapply(fns, function(n) any(calls_of[[n]] %in% vcg_entry), logical(1))
repeat {
  grew <- vapply(fns, function(n) {
    reaches[[n]] || any(reaches[intersect(calls_of[[n]], fns)])
  }, logical(1))
  if (identical(grew, reaches)) break
  reaches <- grew
}

need     <- sort(intersect(getNamespaceExports("ravetools"), fns[reaches]))
ran      <- ls(succeeded)
tried    <- ls(attempted)
missing  <- setdiff(need, tried)              # never called at all
failed   <- setdiff(intersect(need, tried), ran)  # called, but always errored

cat("vcglib-reaching exported functions:", length(need), "\n")
cat("exercised successfully            :", length(intersect(need, ran)), "\n")

status <- 0L
if (length(missing)) {
  cat("\n*** COVERAGE GAP -- never called: ***\n")
  cat(paste0("  - ", missing, collapse = "\n"), "\n")
  status <- 2L
}
if (length(failed)) {
  cat("\n*** COVERAGE GAP -- called but every attempt errored: ***\n")
  cat(paste0("  - ", failed, collapse = "\n"), "\n")
  cat("(a dependency may be missing, or the sample arguments are wrong)\n")
  status <- 2L
}
if (status == 0L) cat("all vcglib-reaching exports exercised\n")

section("done")
if (status != 0L) quit(status = status)
