# Synthetic CCEP-like data with a shared canonical shape plus per-trial scaling.
# `delay` shifts the response onset so the pre-response window is pure noise.
crp_sim <- function(n_time = 600L, n_trials = 30L, delay = 0,
                    t0 = -0.05, t1 = 1.2, noise = 0.25, seed = 1L) {
  set.seed(seed)
  tt <- seq(t0, t1, length.out = n_time)
  peak1 <- delay + 0.10
  peak2 <- delay + 0.30
  shape <- exp(-((tt - peak1) / 0.05)^2) - 0.5 * exp(-((tt - peak2) / 0.10)^2)
  shape[tt < delay] <- 0
  V <- outer(shape, runif(n_trials, 0.6, 1.4)) +
    matrix(rnorm(n_time * n_trials, sd = noise), n_time, n_trials)
  list(V = V * 3, time = tt)
}

test_that("crp returns C_full spanning the entire loaded time, matching C", {
  sim <- crp_sim(seed = 11L)
  res <- crp(sim$V, sim$time)

  parms <- res$parameters
  eps <- sqrt(.Machine$double.eps)

  # C_full spans the entire loaded time range (incl. pre-stimulus t < 0)
  expect_equal(length(parms$C_full), length(sim$time))
  expect_equal(length(parms$params_times_full), length(parms$C_full))
  expect_equal(parms$params_times_full, sim$time)
  expect_equal(res$.data$time, sim$time)

  # reported C (detect_onset = FALSE) is the C_full slice over [t_start, tau_R]
  sel <- parms$params_times_full >= res$t_start &
    parms$params_times_full <= res$tau_R
  expect_equal(parms$C, parms$C_full[sel], tolerance = eps)
  expect_equal(parms$params_times, parms$params_times_full[sel])

  # default C equals the forward unit-norm eigenvector
  expect_equal(sqrt(sum(parms$C^2)), 1, tolerance = eps)

  # tR_sample (window index) matches the reported-C length and tau_R
  tr <- length(parms$C)
  expect_identical(res$projections$tR_sample, tr)
  expect_equal(res$tau_R, parms$params_times[tr])
})

test_that("detect_onset leaves the forward loading (alpha, C_full) untouched", {
  sim <- crp_sim(delay = 0.3, seed = 7L)
  res_off <- crp(sim$V, sim$time, detect_onset = FALSE)
  res_on  <- crp(sim$V, sim$time, detect_onset = TRUE)

  # forward quantities are byte-for-byte identical regardless of onset
  expect_equal(res_on$parameters$al, res_off$parameters$al)
  expect_equal(res_on$parameters$al_p, res_off$parameters$al_p)
  expect_equal(res_on$parameters$Vsnr, res_off$parameters$Vsnr)
  expect_equal(res_on$parameters$expl_var, res_off$parameters$expl_var)
  expect_equal(res_on$parameters$C_full, res_off$parameters$C_full)
  expect_equal(res_on$parameters$ep, res_off$parameters$ep)
  expect_equal(res_on$tau_R, res_off$tau_R)
  expect_identical(res_on$bad_trials, res_off$bad_trials)
})

test_that("detect_onset trims the reported C to [tau_onset, tau_R]", {
  sim <- crp_sim(delay = 0.3, seed = 7L)
  res <- crp(sim$V, sim$time, detect_onset = TRUE)

  expect_true(is.finite(res$tau_onset))
  expect_gte(res$tau_onset, res$t_start)
  expect_lte(res$tau_onset, res$tau_R)
  expect_true(res$tau_onset_lower <= res$tau_onset)
  expect_true(res$tau_onset_upper >= res$tau_onset)

  parms <- res$parameters
  pt <- parms$params_times
  expect_equal(pt[1L], res$tau_onset)
  expect_equal(pt[length(pt)], res$tau_R)
  expect_equal(length(parms$C), length(pt))

  # the trimmed C is exactly the matching slice of C_full
  full_times <- parms$params_times_full
  sel <- full_times >= res$tau_onset & full_times <= res$tau_R
  expect_equal(parms$C, parms$C_full[sel], tolerance = sqrt(.Machine$double.eps))
})

test_that("onset scan recovers a delayed response latency", {
  delay <- 0.30
  sim <- crp_sim(delay = delay, noise = 0.2, seed = 3L)
  res <- crp(sim$V, sim$time, t_start = 0.02, t_end = 1.1,
             detect_onset = TRUE)

  # estimated onset should land near the true delay
  expect_lt(abs(res$tau_onset - delay), 0.1)
  # and well before the response duration tau_R
  expect_lt(res$tau_onset, res$tau_R)
})

test_that("onset can be detected before t_start", {
  delay <- 0.05
  sim <- crp_sim(delay = delay, noise = 0.2, seed = 3L)
  # t_start deliberately set above the true onset
  res <- crp(sim$V, sim$time, t_start = 0.1, t_end = 1.1,
             detect_onset = TRUE)

  expect_true(is.finite(res$tau_onset))
  expect_lt(res$tau_onset, res$t_start)        # earlier than the window start
  expect_lt(abs(res$tau_onset - delay), 0.1)   # near the true latency

  # reported C is sliced from C_full and starts at tau_onset (< t_start)
  parms <- res$parameters
  sel <- parms$params_times_full >= res$tau_onset &
    parms$params_times_full <= res$tau_R
  expect_equal(parms$C, parms$C_full[sel], tolerance = sqrt(.Machine$double.eps))
  expect_lt(parms$params_times[1L], res$t_start)
})

test_that("onset_search_start floors the reverse scan and validates", {
  sim <- crp_sim(delay = 0.05, noise = 0.2, seed = 3L)
  res <- crp(sim$V, sim$time, t_start = 0.1, t_end = 1.1,
             detect_onset = TRUE, onset_search_start = 0.08)
  expect_gte(res$tau_onset, 0.08)

  expect_error(
    crp(sim$V, sim$time, detect_onset = TRUE, onset_search_start = c(0, 1)),
    "onset_search_start"
  )
  expect_error(
    crp(sim$V, sim$time, detect_onset = TRUE, onset_search_start = NA),
    "onset_search_start"
  )
})

test_that("crp drops non-finite time points (NA / NaN / Inf)", {
  sim <- crp_sim(seed = 8L)
  V <- sim$V
  bad <- c(5L, 17L, 100L)
  V[bad[1L], 1L] <- NA
  V[bad[2L], 3L] <- NaN
  V[bad[3L], 2L] <- Inf

  res <- crp(V, sim$time)
  expect_equal(nrow(res$.data$V), length(sim$time) - length(bad))
  expect_equal(length(res$.data$time), length(sim$time) - length(bad))
  expect_false(any(sim$time[bad] %in% res$.data$time))
  expect_equal(length(res$parameters$C_full), length(sim$time) - length(bad))
})

test_that("crp clamps t_start / t_end instead of erroring", {
  sim <- crp_sim(seed = 8L)

  # t_end beyond the last sample -> clamped to max(time), no error
  res <- crp(sim$V, sim$time, t_end = 99)
  expect_equal(res$t_end, max(sim$time))

  # t_start below the first sample -> clamped to min(time)
  res2 <- crp(sim$V, sim$time, t_start = -10, t_end = 1)
  expect_equal(res2$t_start, min(sim$time))

  # impossible window still errors
  expect_error(crp(sim$V, sim$time, t_start = 0.9, t_end = 0.5),
               "t_start.*less than.*t_end")
})

test_that("crp coerces detect_onset silently (isTRUE(as.logical(...)))", {
  sim <- crp_sim(seed = 5L)
  # NA, multi-value, and non-logical all coerce to FALSE — no error
  expect_error(crp(sim$V, sim$time, detect_onset = NA), NA)
  expect_error(crp(sim$V, sim$time, detect_onset = c(TRUE, FALSE)), NA)
  expect_error(crp(sim$V, sim$time, detect_onset = "yes"), NA)
  # "yes" coerces to NA via as.logical, then isTRUE -> FALSE
  r <- crp(sim$V, sim$time, detect_onset = "yes")
  expect_true(is.na(r$tau_onset))
})

test_that("plot.ravetools_crp runs with and without onset", {
  sim <- crp_sim(delay = 0.3, seed = 9L)
  res_off <- crp(sim$V, sim$time, detect_onset = FALSE)
  res_on  <- crp(sim$V, sim$time, detect_onset = TRUE)

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  expect_error(plot(res_off), NA)
  expect_error(plot(res_on), NA)
})
