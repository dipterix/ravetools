# ANTs / ITK-compatible transform I/O

rand_affine <- function(seed = 1) {
  set.seed(seed)
  # a proper rotation (QR of a random matrix) plus an anisotropic scale and shift
  q <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  if (det(q) < 0)
    q[, 1] <- -q[, 1]
  a <- q %*% diag(c(1.05, 0.97, 1.1))
  m <- diag(4)
  m[1:3, 1:3] <- a
  m[1:3, 4] <- c(12.5, -7.25, 3.0)
  m
}

call_rpyANTs <- function(.f, ...) {
  if (system.file(package = "rpyANTs") == "") {
    stop("`rpyANTs` is not installed")
  }
  rpyANTs <- asNamespace("rpyANTs")
  f <- rpyANTs[[.f]]
  f(...)
}

test_that("ras_lps_conjugate is an involution with the expected sign pattern", {
  m <- rand_affine(7)
  expect_equal(ras_lps_conjugate(ras_lps_conjugate(m)), m, tolerance = 1e-12)
  l <- ras_lps_conjugate(m)
  # F = diag(-1,-1,1,1): negate M[1,3], M[1,4], M[2,3], M[2,4], M[3,1], M[3,2]
  expect_equal(l[1, 3], -m[1, 3])
  expect_equal(l[1, 4], -m[1, 4])
  expect_equal(l[2, 3], -m[2, 3])
  expect_equal(l[2, 4], -m[2, 4])
  expect_equal(l[3, 1], -m[3, 1])
  expect_equal(l[3, 2], -m[3, 2])
  # unaffected entries
  expect_equal(l[1, 1], m[1, 1])
  expect_equal(l[3, 3], m[3, 3])
  # 3x4 input is accepted and squared up
  expect_equal(ras_lps_conjugate(m[1:3, ]), l)
})

test_that("affine .mat round-trips at full double precision", {
  m <- rand_affine(2)
  f <- tempfile(fileext = ".mat")
  write_ants_transform(m, f)
  # ITK level-4 .mat: no 128-byte header, two double variables
  expect_equal(unname(file.size(f)), 193)
  expect_equal(read_ants_transform(f), m, tolerance = 1e-12)
})

test_that("warp field round-trips through an ANTs 5-D NIfTI", {
  skip_if_not_installed("freesurferformats")
  set.seed(3)
  nd <- c(7L, 6L, 5L)
  field <- array(rnorm(prod(nd) * 3, sd = 2), c(nd, 3))
  v2r <- matrix(c(0.9, 0, 0, -3, 0, 1.1, 0, -2, 0, 0, 1.2, -1, 0, 0, 0, 1), 4, 4, byrow = TRUE)

  f <- tempfile(fileext = ".nii.gz")
  write_ants_warp(field, f, vox2ras = v2r)
  g <- read_ants_warp(f)
  expect_equal(dim(g), c(nd, 3L))
  expect_equal(as.vector(g), as.vector(field), tolerance = 1e-5)         # float32
  expect_equal(attr(g, "vox2ras"), v2r, tolerance = 1e-5)

  # on disk: 5-D vector image with the vector intent code and our sform
  h <- freesurferformats::read.nifti1.header(f)
  expect_equal(h$dim[1], 5L)
  expect_equal(h$dim[c(2, 3, 4, 5, 6)], c(nd, 1L, 3L))
  expect_equal(h$intent_code, 1007L)
  expect_equal(h$sform_code, 1L)
  expect_equal(h$srow_x, v2r[1, ], tolerance = 1e-5)

  # vox2ras can also travel as an attribute on the field
  attr(field, "vox2ras") <- v2r
  f2 <- tempfile(fileext = ".nii.gz")
  write_ants_warp(field, f2)
  expect_equal(attr(read_ants_warp(f2), "vox2ras"), v2r, tolerance = 1e-5)
})

test_that("the affine can be hidden in / recovered from the warp header (opt-in)", {
  skip_if_not_installed("freesurferformats")
  set.seed(4)
  nd <- c(6L, 5L, 4L)
  field <- array(rnorm(prod(nd) * 3), c(nd, 3))
  v2r <- diag(4)
  v2r[1:3, 4] <- c(-3, -2.5, -2)
  m <- rand_affine(5)

  f <- tempfile(fileext = ".nii.gz")
  write_ants_warp(field,
                  f,
                  vox2ras = v2r,
                  affine = m,
                  direction = "inverse")

  # off by default
  g0 <- read_ants_warp(f)
  expect_null(attr(g0, "transform"))

  # opt-in recovery; lossy precision (5 dp rotation, 2 dp translation)
  g <- read_ants_warp(f, recover_affine = TRUE)
  expect_false(is.null(attr(g, "transform")))
  expect_equal(attr(g, "direction"), "inverse")
  expect_equal(attr(g, "transform")[1:3, 1:3], m[1:3, 1:3], tolerance = 1e-4)
  expect_equal(attr(g, "transform")[1:3, 4], m[1:3, 4], tolerance = 1e-2)

  # values that overflow the fixed-width encoding are skipped with a warning
  big <- diag(4)
  big[1, 1] <- 42
  f2 <- tempfile(fileext = ".nii.gz")
  expect_warning(write_ants_warp(field, f2, vox2ras = v2r, affine = big))
  expect_null(attr(read_ants_warp(f2, recover_affine = TRUE), "transform"))
})

test_that("save_registration / load_registration round-trip via the manifest", {
  skip_if_not_installed("freesurferformats")
  set.seed(6)
  nd <- c(8L, 7L, 6L)
  m <- rand_affine(6)
  v2r <- diag(4)
  v2r[1:3, 4] <- c(-4, -3.5, -3)
  fwd <- array(rnorm(prod(nd) * 3), c(nd, 3))
  attr(fwd, "vox2ras") <- v2r
  inv <- array(rnorm(prod(nd) * 3), c(nd, 3))
  attr(inv, "vox2ras") <- v2r
  res <- structure(
    list(
      transform = m,
      forward_field = fwd,
      inverse_field = inv,
      type = "syn",
      metric = "mattes",
      geometry = list(
        source_vox2ras = v2r,
        target_vox2ras = v2r,
        target_dim = nd
      )
    ),
    class = "ravetools_register_volume3d"
  )

  d <- tempfile()
  dir.create(d)
  man <- save_registration(res, d)
  expect_true(file.exists(man))
  expect_setequal(
    list.files(d),
    c(
      "registration.dcf",
      "registration0GenericAffine.mat",
      "registration1Warp.nii.gz",
      "registration1InverseWarp.nii.gz"
    )
  )

  # load from a different working directory: paths resolve relative to manifest
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  lr <- load_registration(man)
  setwd(old)

  expect_s3_class(lr, "ravetools_register_volume3d")
  expect_equal(lr$transform, m, tolerance = 1e-12)            # .mat is full precision
  expect_equal(as.vector(lr$forward_field), as.vector(fwd), tolerance = 1e-5)
  expect_equal(as.vector(lr$inverse_field), as.vector(inv), tolerance = 1e-5)
  expect_equal(lr$geometry$target_dim, nd)
  expect_equal(lr$geometry$source_vox2ras, v2r, tolerance = 1e-5)
  expect_equal(lr$type, "syn")

  # single-file dispatch
  expect_equal(load_registration(file.path(d, "registration0GenericAffine.mat")), m, tolerance = 1e-12)
  g <- load_registration(file.path(d, "registration1InverseWarp.nii.gz"),
                         recover_affine_from_header = TRUE)
  expect_equal(dim(g), c(nd, 3L))
  expect_equal(attr(g, "direction"), "inverse")
})

# ---- cross-validation against a real ANTs install --------------------------
# All rpyANTs calls should be

ants_ready <- function() {
  isTRUE(tryCatch({
    call_rpyANTs("ants_available")
  }, error = function(e) {
    FALSE
  }))
}

test_that("our affine .mat is read correctly by ANTs (and vice versa)", {
  skip_on_cran()
  skip_if_not_installed("rpyANTs")
  skip_if_not(ants_ready())
  ants <- call_rpyANTs("load_ants")
  m <- rand_affine(11)

  # forward: write ours, read with ANTs. as_ANTsTransform(path)[] is the LPS
  # (ITK-native) matrix, so conjugate back to RAS to compare.
  f <- tempfile(fileext = ".mat")
  write_ants_transform(m, f)
  m_lps <- call_rpyANTs("as_ANTsTransform", f)[]
  expect_equal(ras_lps_conjugate(m_lps), m, tolerance = 1e-5)

  # reverse: ANTs writes the LPS matrix, we read it back as RAS.
  tx <- call_rpyANTs("as_ANTsTransform", ras_lps_conjugate(m))
  f2 <- tempfile(fileext = ".mat")
  ants$write_transform(tx, f2)
  expect_equal(read_ants_transform(f2), m, tolerance = 1e-5)
})

test_that("our warp NIfTI is read correctly by ANTs", {
  skip_on_cran()
  skip_if_not_installed("rpyANTs")
  skip_if_not_installed("freesurferformats")
  skip_if_not(ants_ready())
  ants <- call_rpyANTs("load_ants")

  set.seed(12)
  nd <- c(7L, 6L, 5L)
  field <- array(rnorm(prod(nd) * 3, sd = 1.5), c(nd, 3))
  v2r <- matrix(c(0.9, 0, 0, -3, 0, 1.1, 0, -2, 0, 0, 1.2, -1, 0, 0, 0, 1), 4, 4, byrow = TRUE)
  f <- tempfile(fileext = ".nii.gz")
  write_ants_warp(field, f, vox2ras = v2r)

  img <- ants$image_read(f)
  expect_equal(as.integer(reticulate::py_to_r(img$components)), 3L)
  expect_equal(unlist(reticulate::py_to_r(img$shape)), nd, tolerance = 1e-6)
  expect_equal(unlist(reticulate::py_to_r(img$spacing)), c(0.9, 1.1, 1.2), tolerance = 1e-5)

  # ANTs stores vectors in LPS: x,y components negated relative to our RAS field
  arr <- reticulate::py_to_r(img$numpy())
  expect_equal(dim(arr), c(nd, 3L))
  lps <- field
  lps[, , , 1] <- -lps[, , , 1]
  lps[, , , 2] <- -lps[, , , 2]
  expect_equal(as.vector(arr), as.vector(lps), tolerance = 1e-5)
})
