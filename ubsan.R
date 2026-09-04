# Workload for ubsan.sh -- see that script for how the sanitizer build is made.
#
# Everything runs in ONE R session on purpose. vcglib hands out "user bit"
# flags from function-local statics that live for the whole process, so
# defects of that shape (the CountNonManifoldEdgeFF leak fixed in 0.3.1.1)
# only appear after several calls and are invisible to any single-call test.

library(ravetools)

cat("ravetools:", as.character(packageVersion("ravetools")),
    "from", dirname(system.file(package = "ravetools")), "\n")

section <- function(x) cat("\n========== ", x, " ==========\n", sep = "")

# ---- 1. full test suite ------------------------------------------------
section("testthat::test_local()")
pkg <- Sys.getenv("RAVETOOLS_SRC", unset = ".")
try(testthat::test_local(path = pkg, stop_on_failure = FALSE,
                         reporter = "summary"))

# ---- 2. the vignette CRAN re-builds ------------------------------------
section("vignettes/mesh-geometry.Rmd")
vig <- file.path(pkg, "vignettes", "mesh-geometry.Rmd")
if (file.exists(vig) && requireNamespace("rmarkdown", quietly = TRUE)) {
  out <- file.path(tempdir(), "mesh-geometry.html")
  try(rmarkdown::render(vig, output_file = out, quiet = TRUE,
                        envir = new.env(parent = globalenv())))
} else {
  cat("skipped (vignette or rmarkdown unavailable)\n")
}

# ---- 3. repeat-call stress loop ----------------------------------------
# The point of this sweep: exhaust process-lifetime state. Each entry point is
# called many times in the same session, which is what the 0.3.1.1 bug needed
# (it took four calls to corrupt the face bit allocator).
section("repeat-call stress loop")

N <- 25L
mesh <- vcg_sphere()

# a deliberately defective mesh: duplicating a face makes its three edges
# non-manifold, which drives the topology-repair and diagnostic paths
broken <- mesh
broken$it <- cbind(broken$it, broken$it[, 1, drop = FALSE])

quiet <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(force(expr))),
           error = function(e) invisible(NULL))
}

for (i in seq_len(N)) {
  quiet(vcg_mesh_volume(mesh))
  quiet(vcg_mesh_volume(broken))          # throws; exercises the error path
  quiet(vcg_count_edge_defects(mesh))
  quiet(vcg_count_edge_defects(broken))
  quiet(vcg_average_edge_length(mesh))
  quiet(vcg_max_edge_length(mesh))
  quiet(vcg_update_normals(mesh))
  quiet(vcg_fix_defects(broken, verbose = FALSE))
  quiet(vcg_smooth_explicit(mesh))
  quiet(vcg_smooth_implicit(mesh))
  quiet(vcg_subdivision(mesh))
  quiet(vcg_subdivide_max_edge_length(mesh,
          max_edge_len = vcg_max_edge_length(mesh) * 0.6))
  quiet(vcg_uniform_remesh(mesh, voxel_size = 0.3, verbose = FALSE))
  quiet(vcg_subset_vertex(mesh, seq_len(50)))
  quiet(vcg_kdtree_nearest(mesh$vb[1:3, , drop = FALSE],
                           mesh$vb[1:3, 1:10, drop = FALSE], k = 3))
  quiet(vcg_raycaster(mesh,
                      rbind(mesh$vb[1:3, 1:5, drop = FALSE] * 2,
                            matrix(1, 1, 5)),
                      matrix(-1, 3, 5)))
  quiet(vcg_detect_collision(mesh, mesh))
  quiet(mris_curvature(mesh))
  quiet(mris_smooth(mesh, n_iter = 2))
  quiet(mris_remesh(mesh))
  quiet(mris_inflate(mesh, n_iter = 2))
}
cat("completed", N, "iterations\n")

# isosurface works off a volume, not a mesh
vol <- array(0, c(24, 24, 24))
vol[8:16, 8:16, 8:16] <- 1
for (i in seq_len(N)) {
  quiet(vcg_isosurface(vol, threshold_lb = 0.5))
}
cat("completed", N, "isosurface iterations\n")

section("done")
