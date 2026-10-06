# Acceptance tests for segment_volume_tissue_gmm ("Atropos-lite"), contract D.
#
# Phantom: nested spheres on a 64^3 grid, three classes with means 30/80/130
# (class 1 = outer shell, 2 = middle shell, 3 = inner sphere), Gaussian noise.

gmm_phantom <- function(nd = 64L, noise_sd = 12, seed = 1L,
                        means = c(30, 80, 130), radii = c(28, 20, 12)) {
  set.seed(seed)
  ctr <- (nd - 1) / 2
  idx <- arrayInd(seq_len(nd^3), rep(nd, 3)) - 1
  r <- sqrt(rowSums((idx - ctr)^2))
  truth <- integer(nd^3)
  truth[r <= radii[1]] <- 1L
  truth[r <= radii[2]] <- 2L
  truth[r <= radii[3]] <- 3L
  dim(truth) <- rep(nd, 3)
  vol <- array(0, rep(nd, 3))
  inside <- truth > 0
  vol[inside] <- means[truth[inside]] + stats::rnorm(sum(inside), sd = noise_sd)
  list(volume = vol, truth = truth, mask = inside)
}

# Circular separable Gaussian blur through the FFT (base R only)
gmm_blur3d <- function(x, sigma) {
  d <- dim(x)
  g1 <- function(n) {
    i <- c(0:(n %/% 2), -((n - n %/% 2 - 1):1))
    exp(-i^2 / (2 * sigma^2))
  }
  k <- outer(outer(g1(d[1]), g1(d[2])), g1(d[3]))
  k <- k / sum(k)
  re <- Re(stats::fft(stats::fft(x) * stats::fft(k), inverse = TRUE)) / prod(d)
  re[re < 0] <- 0
  re
}

gmm_dice <- function(a, b) {
  2 * sum(a & b) / (sum(a) + sum(b))
}

# First-maximum argmax over a list of equally shaped arrays
gmm_argmax <- function(lst) {
  best <- lst[[1]]
  am <- rep(1L, length(best))
  for (k in seq_along(lst)[-1]) {
    upd <- lst[[k]] > best
    am[upd] <- k
    best[upd] <- lst[[k]][upd]
  }
  am
}

# b[i] = a[i + off] for every voxel i (zero where i + off leaves the grid), so
# `gmm_shift_neighbors(a, off)[i]` reads the neighbor of i at offset `off`;
# with a negative offset it shifts the content of `a` forward along that axis
gmm_shift_neighbors <- function(a, off) {
  d <- dim(a)
  out <- array(0, d)
  src <- lapply(1:3, function(ax) seq_len(d[ax]) + off[ax])
  ok <- lapply(1:3, function(ax) src[[ax]] >= 1 & src[[ax]] <= d[ax])
  out[ok[[1]], ok[[2]], ok[[3]]] <- a[src[[1]][ok[[1]]], src[[2]][ok[[2]]], src[[3]][ok[[3]]]]
  out
}

# Unbiased weighted variance (ITK's weighted covariance used by Atropos)
gmm_wvar <- function(y, w, m, floor) {
  den <- sum(w) - sum(w^2) / sum(w)
  if (!(den > 0)) return(floor)
  max(sum(w * (y - m)^2) / den, floor)
}

# One EM iteration of the Atropos (Socrates) model with the prior
# initialization, written out in R to pin every detail of the C++ E-step and
# M-step: posterior_k ~ s_k^w (N(y | mu_k, sigma_k^2) M_k)^(1 - w) with the
# spatial prior s_k = pi_k p_k / sum_c pi_c p_c, the normalized mean-field MRF
# term M_k and 1e-10 floors on s_k, M_k and the density
gmm_reference_iteration <- function(volume, mask, priors, prior_weight, mrf_beta,
                                    spacing = c(1, 1, 1), eps = 1e-10) {
  K <- length(priors)
  d <- dim(volume)
  idx <- which(mask)
  y <- volume[idx]
  N <- length(idx)
  P <- vapply(priors, function(p) p[idx], numeric(N))
  s <- rowSums(P)
  informative <- s > 0
  Pn <- matrix(0, N, K)
  Pn[informative, ] <- P[informative, , drop = FALSE] / s[informative]
  lab <- max.col(Pn, ties.method = "first")
  lab[!informative] <- NA_integer_
  var_floor <- 1e-6 * mean((y - mean(y))^2)

  # initial parameters: hard labels weighted by the prior value, label fractions
  mu <- vr <- prop <- numeric(K)
  for (k in seq_len(K)) {
    sel <- which(lab == k)
    wk <- Pn[cbind(sel, k)]
    mu[k] <- sum(wk * y[sel]) / sum(wk)
    vr[k] <- gmm_wvar(y[sel], wk, mu[k], var_floor)
    prop[k] <- length(sel) / sum(informative)
  }

  # neighborhood sums with the one-hot initial labels (unlabeled voxels and
  # voxels outside the mask contribute nothing)
  S <- matrix(0, N, K)
  logM <- matrix(0, N, K)
  if (mrf_beta > 0) {
    onehot <- lapply(seq_len(K), function(k) {
      a <- array(0, d)
      a[idx[which(lab == k)]] <- 1
      a
    })
    for (dz in -1:1) for (dy in -1:1) for (dx in -1:1) {
      if (dx == 0 && dy == 0 && dz == 0) next
      dist <- sqrt(sum((c(dx, dy, dz) * spacing)^2))
      for (k in seq_len(K)) {
        S[, k] <- S[, k] + gmm_shift_neighbors(onehot[[k]], c(dx, dy, dz))[idx] / dist
      }
    }
    logM <- mrf_beta * S
    mx <- apply(logM, 1, max)
    logM <- logM - (mx + log(rowSums(exp(logM - mx))))
    logM <- pmax(logM, log(eps))
  }

  # spatial prior re-weighted by the proportions; zero where all priors vanish
  den <- as.vector(Pn %*% prop)
  sp <- Pn * rep(prop, each = N)
  sp[den > 0, ] <- sp[den > 0, , drop = FALSE] / den[den > 0]
  sp[den <= 0, ] <- Pn[den <= 0, , drop = FALSE]
  log_sp <- log(pmax(sp, eps))
  loglik <- vapply(seq_len(K), function(k) {
    -0.5 * log(2 * base::pi * vr[k]) - (y - mu[k])^2 / (2 * vr[k])
  }, numeric(N))
  loglik <- pmax(loglik, log(eps))
  v <- prior_weight * log_sp + (1 - prior_weight) * (loglik + logM)
  v <- v - apply(v, 1, max)
  post <- exp(v) / rowSums(exp(v))

  # M-step
  mu2 <- colSums(post * y) / colSums(post)
  vr2 <- vapply(seq_len(K), function(k) gmm_wvar(y, post[, k], mu2[k], var_floor), 1)
  list(labels = lab, posteriors = post, means = mu2, sds = sqrt(vr2),
       proportions = colMeans(post), trace = mean(apply(post, 1, max)),
       init_means = mu, init_sds = sqrt(vr), init_proportions = prop)
}

# The Atropos stopping rule on a returned trace: no iteration before the last
# may have stopped (its signed change must be >= tolerance; the first never
# stops), and if fewer than `iterations` ran the last change was < tolerance
gmm_expect_stop_rule <- function(trace, iterations, tolerance) {
  n <- length(trace)
  expect_gte(n, 1L)
  expect_lte(n, iterations)
  if (n > 2L) expect_true(all(diff(trace)[seq_len(n - 2L)] >= tolerance))
  if (n < iterations) {
    expect_gte(n, 2L)
    expect_lt(trace[n] - trace[n - 1L], tolerance)
  }
}

# Deterministic 1-D k-means with equally spaced seeds (the ANTs convention),
# returning the labels of `y` by ascending centroid
gmm_reference_kmeans <- function(y, K) {
  ctr <- min(y) + (max(y) - min(y)) * (seq_len(K) - 0.5) / K
  for (it in 1:200) {
    lab <- max.col(-abs(outer(y, ctr, `-`)), ties.method = "first")
    new <- vapply(seq_len(K), function(k) if (any(lab == k)) mean(y[lab == k]) else ctr[k], 1)
    if (identical(new, ctr)) break
    ctr <- new
  }
  ord <- order(ctr)
  match(lab, ord)
}

test_that("D1: k-means initialization with MRF recovers the phantom; MRF helps at high noise", {
  ph <- gmm_phantom(noise_sd = 12, seed = 1L)
  res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)

  expect_identical(dim(res$segmentation), dim(ph$volume))
  expect_true(is.integer(res$segmentation))
  expect_true(all(res$segmentation[!ph$mask] == 0L))
  expect_true(all(res$segmentation[ph$mask] %in% 1:3))

  # Classes are ordered by ascending mean
  expect_length(res$means, 3L)
  expect_true(all(diff(res$means) > 0))

  for (k in 1:3) {
    expect_gte(gmm_dice(res$segmentation == k, ph$truth == k), 0.95)
  }

  # At noise SD 20 the MRF (beta = 0.2) must beat the plain mixture (beta = 0)
  ph20 <- gmm_phantom(noise_sd = 20, seed = 2L)
  res_mrf <- segment_volume_tissue_gmm(ph20$volume, mask = ph20$mask, mrf_beta = 0.2)
  res_plain <- segment_volume_tissue_gmm(ph20$volume, mask = ph20$mask, mrf_beta = 0)
  dice_mrf <- vapply(1:3, function(k) gmm_dice(res_mrf$segmentation == k, ph20$truth == k), 1)
  dice_plain <- vapply(1:3, function(k) gmm_dice(res_plain$segmentation == k, ph20$truth == k), 1)
  expect_gte(mean(dice_mrf) - mean(dice_plain), 0.03)
})

test_that("D2: priors (blurred truth, shuffled) with prior_weight = 0 keep the prior order", {
  ph <- gmm_phantom(noise_sd = 12, seed = 3L)
  perm <- c(3L, 1L, 2L)
  priors <- lapply(perm, function(k) gmm_blur3d(ph$truth == k, sigma = 2))

  res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                   prior_weight = 0)
  expect_length(res$posteriors, 3L)
  true_means <- c(30, 80, 130)
  for (j in 1:3) {
    expect_gte(gmm_dice(res$segmentation == j, ph$truth == perm[j]), 0.95)
  }
  # Output class j follows prior j (shuffled order), so the means are permuted
  expect_equal(order(res$means), order(perm))
  expect_lt(max(abs(res$means - true_means[perm]) / true_means[perm]), 0.05)
})

test_that("D3: prior_weight = 1 reproduces the Atropos spatial prior argmax(pi * prior)", {
  ph <- gmm_phantom(noise_sd = 12, seed = 4L)
  set.seed(5L)
  flat <- array(stats::rnorm(length(ph$volume), mean = 100, sd = 10), dim(ph$volume))
  priors <- lapply(1:3, function(k) gmm_blur3d(ph$truth == k, sigma = 2))

  # after one iteration the proportions are the initial label fractions
  am0 <- gmm_argmax(priors)
  pi0 <- tabulate(am0[ph$mask], nbins = 3L) / sum(ph$mask)
  res1 <- segment_volume_tissue_gmm(flat, mask = ph$mask, priors = priors,
                                    prior_weight = 1, mrf_beta = 0, iterations = 1L)
  am_pi0 <- gmm_argmax(lapply(1:3, function(k) priors[[k]] * pi0[k]))
  expect_gte(mean(res1$segmentation[ph$mask] == am_pi0[ph$mask]), 0.999)

  # the full run against the returned (converged) proportions
  res <- segment_volume_tissue_gmm(flat, mask = ph$mask, priors = priors,
                                   prior_weight = 1, mrf_beta = 0)
  am_pi <- gmm_argmax(lapply(1:3, function(k) priors[[k]] * res$proportions[k]))
  expect_gte(mean(res$segmentation[ph$mask] == am_pi[ph$mask]), 0.95)

  # at w = 1 neither the image nor the MRF enters the posterior
  res_img <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                       prior_weight = 1, mrf_beta = 0.2)
  expect_identical(res_img$segmentation, res$segmentation)
  expect_equal(res_img$posteriors, res$posteriors, tolerance = 1e-12)
  expect_equal(res_img$proportions, res$proportions, tolerance = 1e-12)
})

test_that("D4: posteriors are proper and the segmentation is their argmax", {
  ph <- gmm_phantom(noise_sd = 12, seed = 6L)
  res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)

  expect_length(res$posteriors, 3L)
  total <- Reduce(`+`, res$posteriors)
  expect_lt(max(abs(total[ph$mask] - 1)), 1e-8)
  for (k in 1:3) {
    p <- res$posteriors[[k]]
    expect_identical(dim(p), dim(ph$volume))
    expect_true(all(p >= 0 & p <= 1))
    expect_true(all(p[!ph$mask] == 0))
  }
  am <- gmm_argmax(res$posteriors)
  expect_true(all(res$segmentation[ph$mask] == am[ph$mask]))

  # proportions are the mean posteriors and sum to one
  expect_equal(sum(res$proportions), 1, tolerance = 1e-10)
  expect_equal(res$proportions,
               vapply(res$posteriors, function(p) mean(p[ph$mask]), 1),
               tolerance = 1e-10)
})

test_that("D5: class means within 5% and SDs within 20% of the truth at noise SD 12", {
  ph <- gmm_phantom(noise_sd = 12, seed = 7L)
  res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)
  true_means <- c(30, 80, 130)
  expect_lt(max(abs(res$means - true_means) / true_means), 0.05)
  expect_lt(max(abs(res$sds - 12) / 12), 0.20)
})

test_that("D6: intensities outside the mask have no effect", {
  ph <- gmm_phantom(noise_sd = 12, seed = 8L)
  res_a <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)

  vol_b <- ph$volume
  set.seed(9L)
  outside <- which(!ph$mask)
  vol_b[outside] <- stats::rnorm(length(outside), mean = 5000, sd = 1000)
  vol_b[outside[1:100]] <- NA_real_
  vol_b[outside[101:200]] <- Inf
  vol_b[outside[201:300]] <- -Inf
  res_b <- segment_volume_tissue_gmm(vol_b, mask = ph$mask)

  expect_identical(res_a$segmentation, res_b$segmentation)
  expect_identical(res_a$posteriors, res_b$posteriors)
  expect_identical(res_a$means, res_b$means)
  expect_identical(res_a$sds, res_b$sds)
  expect_identical(res_a$proportions, res_b$proportions)
  expect_identical(res_a$trace, res_b$trace)

  # priors outside the mask are ignored as well
  priors <- lapply(1:3, function(k) gmm_blur3d(ph$truth == k, sigma = 2))
  res_p1 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                      prior_weight = 0.25)
  priors2 <- lapply(priors, function(p) {
    p[!ph$mask] <- 7
    p
  })
  res_p2 <- segment_volume_tissue_gmm(vol_b, mask = ph$mask, priors = priors2,
                                      prior_weight = 0.25)
  expect_identical(res_p1$segmentation, res_p2$segmentation)
  expect_identical(res_p1$posteriors, res_p2$posteriors)
})

test_that("D7: results are deterministic across thread counts", {
  ph <- gmm_phantom(noise_sd = 12, seed = 10L)
  priors <- lapply(1:3, function(k) gmm_blur3d(ph$truth == k, sigma = 2))

  # restore whatever thread setting was active before (an unset variable is "auto")
  old_threads <- Sys.getenv("RAVETOOLS_NUM_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old_threads)) ravetools_threads() else ravetools_threads(n_threads = as.integer(old_threads))
  }, add = TRUE)
  ravetools_threads(n_threads = 1L)
  res_1 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)
  res_1p <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                      prior_weight = 0.25)
  # CRAN allows at most 2 cores; use more only off CRAN
  n_threads <- if (is_not_cran(if_interactive = FALSE)) max(2L, min(8L, detect_threads())) else 2L
  ravetools_threads(n_threads = n_threads)
  res_n <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)
  res_np <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                      prior_weight = 0.25)

  expect_identical(res_1$segmentation, res_n$segmentation)
  expect_lt(max(abs(res_1$means - res_n$means)), 1e-10)
  expect_lt(max(abs(res_1$sds - res_n$sds)), 1e-10)
  for (k in 1:3) {
    expect_lt(max(abs(res_1$posteriors[[k]] - res_n$posteriors[[k]])), 1e-10)
  }
  expect_identical(res_1p$segmentation, res_np$segmentation)
  expect_lt(max(abs(res_1p$means - res_np$means)), 1e-10)
  for (k in 1:3) {
    expect_lt(max(abs(res_1p$posteriors[[k]] - res_np$posteriors[[k]])), 1e-10)
  }
})

test_that("D8: invalid priors, dimensions and parameters error", {
  ph <- gmm_phantom(nd = 24L, noise_sd = 12, seed = 11L, radii = c(10, 7, 4))
  priors <- lapply(1:3, function(k) gmm_blur3d(ph$truth == k, sigma = 2))

  # explicit n_classes conflicting with the number of priors
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                         n_classes = 4L))
  # the default n_classes does not conflict (K comes from the priors)
  expect_silent(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                          iterations = 1L))
  # a prior on the wrong grid
  bad <- priors
  bad[[2]] <- bad[[2]][1:20, , ]
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = bad))
  # negative or missing prior values
  bad <- priors
  bad[[1]][ph$mask][1] <- -1
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = bad))
  bad <- priors
  bad[[3]][ph$mask][1] <- NA
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = bad))
  # too few priors / not a list
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors[1]))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors[[1]]))

  # mask on the wrong grid, non-3D volume, empty mask
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask[1:20, , ]))
  expect_error(segment_volume_tissue_gmm(ph$volume[, , 1], mask = ph$mask[, , 1]))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = array(FALSE, dim(ph$volume))))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, n_classes = 1L))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, n_classes = 2.5))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, iterations = 1.5))
  # n_classes = NULL counts as not given: fine with priors, an error without
  expect_silent(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                          n_classes = NULL, iterations = 1L))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, n_classes = NULL),
               "n_classes")

  # parameter ranges
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, prior_weight = 1.5))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = -1))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = c(1L, 1L)))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = -1L))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, iterations = 0L))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, tolerance = -1))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, vox2ras = diag(3)))

  # non-integer radii are rejected instead of silently truncated
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = 1.9),
               "mrf_radius")
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = 0.9),
               "mrf_radius")
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = c(1, 1.5, 1)),
               "mrf_radius")
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = c(1L, NA, 1L)),
               "mrf_radius")
  # non-finite or out-of-range iterations fail in R, without a coercion warning
  expect_silent(expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, iterations = Inf),
                             "segment_volume_tissue_gmm.*iterations"))
  expect_silent(expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, iterations = 1e10),
                             "segment_volume_tissue_gmm.*iterations"))
  # the same for n_classes (beyond the integer range, non-finite, NA)
  expect_silent(expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, n_classes = 2^31),
                             "segment_volume_tissue_gmm.*n_classes"))
  expect_silent(expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, n_classes = 1e10),
                             "segment_volume_tissue_gmm.*n_classes"))
  expect_silent(expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, n_classes = Inf),
                             "segment_volume_tissue_gmm.*n_classes"))
  expect_silent(expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, n_classes = NA_real_),
                             "segment_volume_tissue_gmm.*n_classes"))
  # a malformed n_classes given with priors is reported as malformed, and only
  # a valid value that differs from the number of priors as a conflict
  e <- tryCatch(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                          n_classes = "3"), error = identity)
  expect_match(conditionMessage(e), "n_classes.*single integer")
  expect_false(grepl("conflicts", conditionMessage(e), fixed = TRUE))
  e <- tryCatch(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                          n_classes = 2.5), error = identity)
  expect_match(conditionMessage(e), "n_classes.*single integer")
  e <- tryCatch(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                          n_classes = 4), error = identity)
  expect_match(conditionMessage(e), "n_classes.*conflicts with the number of priors")
  # an MRF factor so large that beta times the neighborhood weights would
  # overflow is rejected, while a merely huge one still gives proper posteriors
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = 1e308),
               "mrf_beta")
  res_hard <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = 1e300,
                                        iterations = 2L)
  total <- Reduce(`+`, res_hard$posteriors)
  expect_true(all(is.finite(total[ph$mask])))
  expect_lt(max(abs(total[ph$mask] - 1)), 1e-8)
  expect_true(all(res_hard$segmentation[ph$mask] %in% 1:3))
  # numeric arguments must be numeric (no silent coercion of strings or logicals)
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                         prior_weight = "0.5"), "prior_weight")
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = "1"),
               "mrf_radius")
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = "0.2"),
               "mrf_beta")
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, tolerance = "0"),
               "tolerance")
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, iterations = TRUE),
               "iterations")

  # constant intensities inside the mask cannot be clustered
  flat <- ph$volume
  flat[ph$mask] <- 42
  expect_error(segment_volume_tissue_gmm(flat, mask = ph$mask))
})

test_that("mask = NULL uses every finite voxel; vox2ras is propagated; radius 0 equals beta 0", {
  ph <- gmm_phantom(nd = 32L, noise_sd = 12, seed = 12L, radii = c(14, 10, 6))
  vol_na <- ph$volume
  vol_na[!ph$mask] <- NA_real_

  res_mask <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)
  res_null <- segment_volume_tissue_gmm(vol_na)
  expect_identical(res_mask$segmentation, res_null$segmentation)
  expect_identical(res_mask$posteriors, res_null$posteriors)
  expect_null(attr(res_mask$segmentation, "vox2ras"))

  # vox2ras from the argument and from the attribute, propagated to the outputs
  v2r <- diag(c(0.8, 1, 1.5, 1))
  v2r[1:3, 4] <- c(-10, -20, -30)
  res_v <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, vox2ras = v2r)
  expect_equal(attr(res_v$segmentation, "vox2ras"), v2r)
  expect_equal(attr(res_v$posteriors[[2]], "vox2ras"), v2r)
  vol_attr <- ph$volume
  attr(vol_attr, "vox2ras") <- v2r
  res_a <- segment_volume_tissue_gmm(vol_attr, mask = ph$mask)
  expect_identical(res_a$segmentation, res_v$segmentation)
  expect_identical(res_a$posteriors, res_v$posteriors)
  # anisotropic spacing changes the MRF weights, so the posteriors differ
  expect_false(isTRUE(all.equal(res_v$posteriors[[1]], res_mask$posteriors[[1]])))
  for (k in 1:3) {
    expect_gte(gmm_dice(res_v$segmentation == k, ph$truth == k), 0.95)
  }

  # an empty neighborhood is the same model as a zero smoothing factor
  res_r0 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = 0L)
  res_b0 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = 0)
  expect_identical(res_r0$segmentation, res_b0$segmentation)
  expect_equal(res_r0$posteriors, res_b0$posteriors, tolerance = 1e-12)

  # a per-axis radius is accepted; integer and logical inputs are accepted
  res_r3 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = c(1L, 1L, 0L))
  expect_true(all(res_r3$segmentation[ph$mask] %in% 1:3))
  vol_int <- ph$volume
  storage.mode(vol_int) <- "integer"
  res_i <- segment_volume_tissue_gmm(vol_int, mask = ph$mask * 1L, iterations = 2L)
  expect_true(all(res_i$segmentation[ph$mask] %in% 1:3))

  # trailing singleton dimensions are dropped; the outputs are always 3D
  vol4 <- array(ph$volume, c(dim(ph$volume), 1L))
  mask4 <- array(ph$mask, c(dim(ph$mask), 1L, 1L))
  res_4 <- segment_volume_tissue_gmm(vol4, mask = mask4)
  expect_identical(dim(res_4$segmentation), dim(ph$volume))
  expect_identical(dim(res_4$posteriors[[1]]), dim(ph$volume))
  expect_identical(res_4$segmentation, res_mask$segmentation)
  expect_error(segment_volume_tissue_gmm(array(ph$volume, c(dim(ph$volume) / c(1, 1, 2), 2L))))

  # the trace has one entry per iteration run (at most `iterations`, fewer
  # when the Atropos stopping rule fires) and a large tolerance stops early
  gmm_expect_stop_rule(res_mask$trace, iterations = 5L, tolerance = 0)
  expect_true(all(res_mask$trace > 0 & res_mask$trace <= 1))
  res_tol <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, tolerance = 1)
  expect_length(res_tol$trace, 2L)
})

test_that("tolerance follows the Atropos stopping rule: tolerance = 0 stops at the first decrease", {
  ph <- gmm_phantom(nd = 32L, noise_sd = 12, seed = 1L, radii = c(14, 10, 6))

  # without the MRF the mean maximum posterior decreases from the first
  # iteration on this phantom, so the default tolerance = 0 stops after the
  # second iteration ...
  res0 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = 0)
  expect_length(res0$trace, 2L)
  expect_lt(res0$trace[2], res0$trace[1])
  gmm_expect_stop_rule(res0$trace, iterations = 5L, tolerance = 0)
  # ... and returns the state of the iteration in which the decrease happened
  # (no rollback), which is not the state after a single iteration
  res2 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = 0, iterations = 2L)
  expect_identical(res0, res2)
  res1 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = 0, iterations = 1L)
  expect_length(res1$trace, 1L)
  expect_equal(res1$trace, res0$trace[1])
  expect_false(isTRUE(all.equal(res1$means, res0$means)))
  expect_false(identical(res1$segmentation, res0$segmentation))

  # with the MRF the measure increases at every iteration here, so all run
  res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)
  expect_length(res$trace, 5L)
  expect_true(all(diff(res$trace) > 0))
  gmm_expect_stop_rule(res$trace, iterations = 5L, tolerance = 0)
  # a positive tolerance stops as soon as the increase falls below it: the
  # increases on this phantom are about 9e-4, 2.5e-4, 7e-5 and 2e-5
  res_tol <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, tolerance = 1e-4)
  expect_length(res_tol$trace, 4L)
  expect_identical(res_tol$trace, res$trace[1:4])
  gmm_expect_stop_rule(res_tol$trace, iterations = 5L, tolerance = 1e-4)
  res_it4 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, iterations = 4L)
  expect_identical(res_tol, res_it4)

  # the first iteration never stops, whatever the tolerance
  res_big <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, tolerance = 1)
  expect_length(res_big$trace, 2L)
  res_big1 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, tolerance = 1, iterations = 1L)
  expect_length(res_big1$trace, 1L)
  # `iterations` is only an upper bound: a huge value with a tolerance stops
  # right away and gives the same result
  res_huge <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask,
                                        iterations = .Machine$integer.max, tolerance = 1)
  expect_identical(res_huge, res_big)

  # a negative tolerance is invalid input
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, tolerance = -1e-12),
               "tolerance")
})

test_that("errors from the native kernel are attributed to segment_volume_tissue_gmm", {
  ph <- gmm_phantom(nd = 16L, noise_sd = 12, seed = 22L, radii = c(7, 5, 3))

  flat <- ph$volume
  flat[ph$mask] <- 42
  e <- tryCatch(segment_volume_tissue_gmm(flat, mask = ph$mask), error = identity)
  expect_s3_class(e, "error")
  expect_match(conditionMessage(e),
               "^`segment_volume_tissue_gmm`: Intensities inside the mask are constant")
  expect_identical(conditionCall(e)[[1]], as.name("segment_volume_tissue_gmm"))

  vol_hot <- ph$volume
  vol_hot[which(ph$mask)[1]] <- 1e6
  e <- tryCatch(segment_volume_tissue_gmm(vol_hot, mask = ph$mask), error = identity)
  expect_match(conditionMessage(e), "^`segment_volume_tissue_gmm`: Class .* extreme values")
  expect_identical(conditionCall(e)[[1]], as.name("segment_volume_tissue_gmm"))

  # R-level validation errors look the same
  e <- tryCatch(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = -1),
                error = identity)
  expect_match(conditionMessage(e), "^`segment_volume_tissue_gmm`: `mrf_beta`")
  expect_identical(conditionCall(e)[[1]], as.name("segment_volume_tissue_gmm"))
})

test_that("a large common intensity offset does not degrade the initialization", {
  ph <- gmm_phantom(nd = 32L, noise_sd = 12, seed = 21L, radii = c(14, 10, 6))
  res0 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask)
  # at offset 1e15 the shifted intensities are still exact to 0.125, but the
  # raw running sums would exceed 2^53 times the signal and discard it
  off <- 1e15
  res_off <- segment_volume_tissue_gmm(ph$volume + off, mask = ph$mask)
  expect_lt(max(abs((res_off$means - off) - res0$means)), 0.5)
  expect_lt(max(abs(res_off$sds / res0$sds - 1)), 0.02)
  expect_gte(mean(res_off$segmentation[ph$mask] == res0$segmentation[ph$mask]), 0.99)
  expect_lt(max(abs(res_off$trace - res0$trace)), 1e-3)

  # an intensity range whose squared deviations overflow is an error rather
  # than NaN posteriors
  expect_error(segment_volume_tissue_gmm(ph$volume * 1e160, mask = ph$mask), "rescale")
})

test_that("priors that carry no information inside the mask error", {
  ph <- gmm_phantom(nd = 24L, noise_sd = 12, seed = 14L, radii = c(10, 7, 4))

  # all priors zero everywhere
  pz <- lapply(1:3, function(k) array(0, dim(ph$volume)))
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = pz),
               "zero at every voxel")
  # integer-truncated probability maps are all zero as well
  pint <- lapply(1:3, function(k) {
    p <- (ph$truth == k) * 0.9
    storage.mode(p) <- "integer"
    p
  })
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = pint),
               "zero at every voxel")
  # priors that are positive only outside the mask carry no information either
  pout <- lapply(1:3, function(k) {
    p <- gmm_blur3d(ph$truth == k, sigma = 2)
    p[ph$mask] <- 0
    p
  })
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = pout),
               "zero at every voxel")
  # a single class whose prior vanishes inside the mask is reported by class
  pone <- lapply(1:3, function(k) gmm_blur3d(ph$truth == k, sigma = 2))
  pone[[2]][] <- 0
  expect_error(segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = pone),
               "prior 2")

  # isolated zero voxels fall back to a uniform prior and still segment
  piso <- lapply(1:3, function(k) gmm_blur3d(ph$truth == k, sigma = 2))
  zero_at <- which(ph$mask)[1:50]
  for (k in 1:3) piso[[k]][zero_at] <- 0
  res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = piso,
                                   prior_weight = 0.25)
  expect_true(all(res$segmentation[ph$mask] %in% 1:3))
  for (k in 1:3) {
    expect_gte(gmm_dice(res$segmentation == k, ph$truth == k), 0.85)
  }
})

test_that("k-means initialization reports extreme outlier intensities clearly", {
  ph <- gmm_phantom(nd = 24L, noise_sd = 12, seed = 15L, radii = c(10, 7, 4))
  vol_hot <- ph$volume
  vol_hot[which(ph$mask)[1]] <- 1e6   # one hot voxel
  expect_error(segment_volume_tissue_gmm(vol_hot, mask = ph$mask), "extreme values")

  # clipping the outlier restores a normal run
  vol_clip <- vol_hot
  vol_clip[vol_clip > 200] <- 200
  res <- segment_volume_tissue_gmm(vol_clip, mask = ph$mask)
  for (k in 1:3) {
    expect_gte(gmm_dice(res$segmentation == k, ph$truth == k), 0.9)
  }
})

test_that("mrf_radius beyond the grid extent is clamped without changing results", {
  ph <- gmm_phantom(nd = 12L, noise_sd = 12, seed = 16L, radii = c(5, 3.5, 2))
  res_full <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = 11L,
                                        iterations = 2L)
  res_big <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = 20L,
                                       iterations = 2L)
  res_mix <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask,
                                       mrf_radius = c(11L, 40L, 1000L), iterations = 2L)
  res_huge <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = 1e12,
                                        iterations = 2L)
  expect_identical(res_big, res_full)
  expect_identical(res_mix, res_full)
  expect_identical(res_huge, res_full)
  # a radius inside the grid is still a different model
  res_1 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_radius = 1L,
                                     iterations = 2L)
  expect_false(identical(res_1$posteriors, res_full$posteriors))
})

test_that("integer and logical priors through the C++ entry point are safe under gctorture", {
  nd <- 12L
  d <- rep(nd, 3L)
  set.seed(17L)
  truth <- sample.int(3L, nd^3, replace = TRUE)
  vol <- truth * 50 + stats::rnorm(nd^3, sd = 5)
  mask <- rep(TRUE, nd^3)
  pr_dbl <- lapply(1:3, function(k) {
    p <- rep(1, nd^3)
    p[truth == k] <- 1000
    p
  })
  pr_int <- lapply(pr_dbl, function(p) {
    storage.mode(p) <- "integer"
    p
  })
  pr_lgl <- lapply(1:3, function(k) truth == k)
  call_cpp <- function(pr) {
    segment_gmm_mrf_cpp(volume = vol, dims = d, mask = mask, priors = pr, n_classes = 3L,
                        prior_weight = 1, mrf_beta = 0, mrf_radius = c(1L, 1L, 1L),
                        direction = diag(3), iterations = 1L, tolerance = 0, verbose = FALSE)
  }
  ref <- call_cpp(pr_dbl)
  expect_identical(ref$segmentation, truth)
  ref_lgl <- call_cpp(lapply(pr_lgl, as.double))

  # Every allocation triggers a collection: a pointer kept past the scoped
  # handle of a coerced (integer or logical) prior would read freed memory
  gctorture(TRUE)
  on.exit(gctorture(FALSE), add = TRUE)
  res_int <- call_cpp(pr_int)
  res_lgl <- call_cpp(pr_lgl)
  gctorture(FALSE)

  expect_identical(res_int$segmentation, ref$segmentation)
  expect_identical(res_int$posteriors, ref$posteriors)
  expect_identical(res_int$means, ref$means)
  expect_identical(res_lgl$segmentation, ref_lgl$segmentation)
  expect_identical(res_lgl$posteriors, ref_lgl$posteriors)
})

test_that("D10: the pull of partially disagreeing priors increases with prior_weight", {
  ph <- gmm_phantom(noise_sd = 12, seed = 18L)
  # priors shifted by 4 voxels along the first axis disagree with the image
  # along every tissue boundary
  priors <- lapply(1:3, function(k) {
    gmm_shift_neighbors(gmm_blur3d(ph$truth == k, sigma = 2), c(-4L, 0L, 0L))
  })
  am <- gmm_argmax(priors)
  ws <- c(0, 0.25, 0.5, 0.75)
  runs <- lapply(ws, function(w) {
    segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors, prior_weight = w)
  })
  agree_prior <- vapply(runs, function(r) mean(r$segmentation[ph$mask] == am[ph$mask]), 1)
  dice_truth <- vapply(runs, function(r) {
    mean(vapply(1:3, function(k) gmm_dice(r$segmentation == k, ph$truth == k), 1))
  }, 1)

  # agreement with the priors grows with the weight, and materially so
  expect_true(all(diff(agree_prior) > 0))
  expect_gt(agree_prior[4] - agree_prior[1], 0.05)
  # with w = 0 the priors only initialize and the image wins; at w = 0.75 the
  # shifted priors have pulled the segmentation well away from the truth
  expect_gte(dice_truth[1], 0.95)
  expect_lt(dice_truth[4], dice_truth[1] - 0.1)
})

test_that("one EM iteration matches an R reference of the Atropos (Socrates) model", {
  ph <- gmm_phantom(nd = 24L, noise_sd = 12, seed = 19L, radii = c(10, 7, 4))
  perm <- c(2L, 3L, 1L)
  priors <- lapply(perm, function(k) gmm_blur3d(ph$truth == k, sigma = 2))
  # scattered voxels where every prior is zero start unlabeled
  zero_at <- which(ph$mask)[seq(1, sum(ph$mask), by = 97)]
  for (k in 1:3) priors[[k]][zero_at] <- 0
  v2r <- diag(c(1, 1.5, 2, 1))
  spacing <- c(1, 1.5, 2)

  for (w in c(0, 0.5, 1)) {
    for (beta in c(0, 0.2)) {
      tag <- sprintf("w = %g, beta = %g", w, beta)
      ref <- gmm_reference_iteration(ph$volume, ph$mask, priors, w, beta, spacing = spacing)
      res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, priors = priors,
                                       prior_weight = w, mrf_beta = beta,
                                       iterations = 1L, vox2ras = v2r)
      for (k in 1:3) {
        p <- res$posteriors[[k]][ph$mask]
        expect_lt(max(abs(p - ref$posteriors[, k])), 1e-8, label = tag)
        # the floors make every posterior positive; their logs must agree too
        expect_lt(max(abs(log(p) - log(ref$posteriors[, k]))), 1e-6, label = tag)
      }
      expect_equal(res$means, ref$means, tolerance = 1e-8, label = tag)
      expect_equal(res$sds, ref$sds, tolerance = 1e-8, label = tag)
      expect_equal(res$proportions, ref$proportions, tolerance = 1e-8, label = tag)
      expect_equal(res$trace, ref$trace, tolerance = 1e-8, label = tag)
    }
  }
})

test_that("without priors the posterior is the floored density times the MRF term, free of proportions", {
  ph <- gmm_phantom(nd = 24L, noise_sd = 12, seed = 20L, radii = c(10, 7, 4))
  y <- ph$volume[ph$mask]
  lab <- gmm_reference_kmeans(y, 3L)
  var_floor <- 1e-6 * mean((y - mean(y))^2)
  mu <- vapply(1:3, function(k) mean(y[lab == k]), 1)
  vr <- vapply(1:3, function(k) gmm_wvar(y[lab == k], rep(1, sum(lab == k)), mu[k], var_floor), 1)
  loglik <- vapply(1:3, function(k) {
    -0.5 * log(2 * base::pi * vr[k]) - (y - mu[k])^2 / (2 * vr[k])
  }, numeric(length(y)))
  loglik <- pmax(loglik, log(1e-10))
  v <- loglik - apply(loglik, 1, max)
  post <- exp(v) / rowSums(exp(v))

  res <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, mrf_beta = 0, iterations = 1L)
  for (k in 1:3) {
    p <- res$posteriors[[k]][ph$mask]
    expect_lt(max(abs(p - post[, k])), 1e-8)
    expect_lt(max(abs(log(p) - log(post[, k]))), 1e-6)
  }
  # the label fractions are far from uniform, so including them would show
  prop0 <- tabulate(lab, 3L) / length(lab)
  expect_gt(max(prop0) - min(prop0), 0.3)

  # prior_weight has no effect without priors
  res_w0 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, prior_weight = 0)
  res_w9 <- segment_volume_tissue_gmm(ph$volume, mask = ph$mask, prior_weight = 0.9)
  expect_identical(res_w0, res_w9)
})
