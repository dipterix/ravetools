
# Tests for the firls() fast-path fix (fullband check + Cholesky solver).
#
# Two invariants to verify:
#   1. Correctness: new path produces coefficients numerically identical to
#      legacy (within ~1e-12 for well-conditioned cases, noting the two paths
#      are algebraically equivalent when fullband=TRUE).
#   2. fullband detection fires correctly:
#      - fir1-style specs (no gaps)  -> fullband = TRUE  -> fast path
#      - specs with inter-band gaps  -> fullband = FALSE -> matrix path
#
# Helper: normalised frequency gain at f0 in [0, 1]
fir_gain <- function(b, f0) {
  L <- length(b)
  Mod(sum(exp(-1i * 2 * pi * seq(0, L - 1) * (f0 / 2)) * b))
}

# -- 1. fir1-compatible specs: no inter-band gaps ------------------------------

test_that("firls new vs legacy agree for no-gap (fir1-compatible) specs", {

  cases <- list(
    # name,  N,   freq,                    A
    list("low-pass even",     100, c(0, 0.3, 0.3, 1),       c(1, 1, 0, 0)),
    list("low-pass odd",      101, c(0, 0.3, 0.3, 1),       c(1, 1, 0, 0)),
    list("high-pass even",    80,  c(0, 0.6, 0.6, 1),       c(0, 0, 1, 1)),
    list("band-pass even",    120, c(0, 0.2, 0.2, 0.5, 0.5, 1), c(0, 0, 1, 1, 0, 0)),
    list("band-stop even",    120, c(0, 0.2, 0.2, 0.5, 0.5, 1), c(1, 1, 0, 0, 1, 1)),
    list("multi-band even",   60,  c(0, 0.1, 0.1, 0.3, 0.3, 0.6, 0.6, 1),
                                   c(1, 1, 0, 0, 1, 1, 0, 0))
  )

  for (case in cases) {
    label <- case[[1]]
    N <- case[[2]]
    freq <- case[[3]]
    A <- case[[4]]

    b_new    <- firls(N, freq, A)$b
    b_legacy <- firls(N, freq, A, legacy = TRUE)$b

    expect_equal(
      b_new, b_legacy, tolerance = 1e-12,
      label = sprintf("firls no-gap case '%s'", label)
    )
  }
})

# -- 2. fir1() output unchanged ------------------------------------------------
# fir1 calls firls internally; verify that the window-based filter coefficients
# are unchanged across several filter types.

test_that("fir1 coefficients are unchanged by the fast-path fix", {

  for (n in c(20, 50, 101)) {
    b_new    <- fir1(n, 0.3, "low")$b
    b_legacy <- firls(n, c(0, 0.3, 0.3, 1), c(1, 1, 0, 0), legacy = TRUE)$b
    # fir1 applies a hamming window on top, so compare the raw firls output
    # rather than the final windowed b.  Instead just check the gain is the same.
    # (fir1 always applies the hamming window, so we re-derive for legacy too)
    wind <- hamming(n + 1)
    b_legacy_win <- b_legacy * wind
    # Normalise both at DC the same way (scale=TRUE in fir1 = b/sum(b) for LP)
    b_legacy_win <- b_legacy_win / sum(b_legacy_win)

    expect_equal(
      b_new, b_legacy_win, tolerance = 1e-12,
      label = sprintf("fir1 low-pass n=%d", n)
    )
  }
})

# -- 3. Gain correctness after fix --------------------------------------------
# Independent check: filters produced by the new path must satisfy their
# passband / stopband specifications to within reasonable tolerance.

test_that("firls (new path) satisfies passband gain requirements", {

  # Low-pass: gain near 1 in passband, near 0 in stopband
  b <- firls(100, c(0, 0.3, 0.3, 1), c(1, 1, 0, 0))$b
  expect_gt(fir_gain(b, 0),   0.99)   # DC gain close to 1
  expect_gt(fir_gain(b, 0.15), 0.9)   # mid-passband
  expect_lt(fir_gain(b, 0.5),  0.05)  # well into stopband

  # High-pass
  b <- firls(100, c(0, 0.6, 0.6, 1), c(0, 0, 1, 1))$b
  expect_lt(fir_gain(b, 0.1),  0.05)
  expect_gt(fir_gain(b, 0.8),  0.9)
  expect_gt(fir_gain(b, 1.0),  0.99)

  # Band-pass
  b <- firls(120, c(0, 0.2, 0.2, 0.6, 0.6, 1), c(0, 0, 1, 1, 0, 0))$b
  expect_lt(fir_gain(b, 0.05), 0.1)
  expect_gt(fir_gain(b, 0.4),  0.9)
  expect_lt(fir_gain(b, 0.85), 0.1)
})

# -- 4. Specs with inter-band gaps must still give correct results -------------
# When freq has genuine gaps (transition bands), fullband=FALSE -> matrix path.
# Verify the Cholesky solver gives the same answer as QR for these cases.

test_that("firls with inter-band gaps: Cholesky agrees with legacy QR", {

  gap_cases <- list(
    list("gap low-pass",  80,  c(0, 0.25, 0.35, 1),       c(1, 1, 0, 0)),
    list("gap band-pass", 80,  c(0, 0.15, 0.25, 0.5, 0.6, 1), c(0, 0, 1, 1, 0, 0)),
    list("gap notch",     80,  c(0, 0.2, 0.3, 0.5, 0.6, 1),   c(1, 1, 0, 0, 1, 1))
  )

  for (case in gap_cases) {
    label <- case[[1]]
    N <- case[[2]]
    freq <- case[[3]]
    A <- case[[4]]

    b_new    <- firls(N, freq, A)$b
    b_legacy <- firls(N, freq, A, legacy = TRUE)$b

    # Both answers are within the backward-error bound of the ill-conditioned
    # system; 1e-10 is tighter than the actual error floor for small N.
    expect_equal(
      b_new, b_legacy, tolerance = 1e-10,
      label = sprintf("firls gap case '%s'", label)
    )
  }
})

# -- 5. design_filter_fir ("firls" method) coefficients unchanged by fix -------
# design_filter_fir builds the firls freq spec as:
#   w   = c(0, rep(bands, each = 2), 1)       # no gaps between adjacent bands
#   a   = rep(mag2, each = 2)[-c(1, end)]
# then applies a hamming window and scales.  The test reproduces that exact call
# with legacy=TRUE and confirms bit-for-bit agreement (within 1e-12).

test_that("design_filter_fir method='firls' coefficients unchanged by fix", {

  fs <- 500

  reconstruct_legacy <- function(f_new) {
    # Recover bands and mag2 from stored parameters
    freq_stored <- f_new$parameters$frequency   # c(0, bands..., 1)
    amp_stored  <- f_new$parameters$amplitude   # mag2
    bands <- freq_stored[-c(1, length(freq_stored))]   # drop 0 and 1
    mag2  <- amp_stored

    # Reproduce design_filter_fir's internal firls call exactly
    w <- c(0, rep(bands, each = 2), 1)
    a <- rep(mag2, each = 2)
    a <- a[-c(1, length(a))]

    b <- firls(N = f_new$parameters$order, freq = w, A = a, legacy = TRUE)$b

    # Same hamming window as design_filter_fir
    wind <- hamming(length(b))
    b <- b * wind

    # Reproduce design_filter_fir's reference-frequency scaling (matches the
    # outer `if (scale)` block in design_filter_fir, not fir1's scale_filter):
    #   low-pass  -> f0 = 0  (DC,      b / sum(b))
    #   high-pass -> f0 = 1  (Nyquist)
    #   band-pass -> f0 = mean(kaisprm$Wc)
    #                      = mean(c(passband, stopband))  [for the firls method]
    #   band-stop -> f0 = 0
    ftype <- f_new$parameters$type[[2]]
    f0 <- switch(
      ftype,
      "low"  = 0,
      "high" = 1,
      "pass" = mean(c(f_new$parameters$passband, f_new$parameters$stopband)),
      0   # stop / default -> DC
    )
    L <- length(b)
    b / Mod(sum(exp(-1i * 2 * pi * seq(0, L - 1) * (f0 / 2)) * b))
  }

  # Low-pass
  f_lp <- design_filter_fir(
    sample_rate = fs, low_pass_freq = 50,
    filter_order = 100, method = "firls"
  )
  expect_equal(f_lp$b, reconstruct_legacy(f_lp), tolerance = 1e-12,
               label = "design_filter_fir firls low-pass coeff equality")

  # Band-pass
  f_bp <- design_filter_fir(
    sample_rate = fs, high_pass_freq = 80, low_pass_freq = 150,
    filter_order = 100, method = "firls"
  )
  expect_equal(f_bp$b, reconstruct_legacy(f_bp), tolerance = 1e-12,
               label = "design_filter_fir firls band-pass coeff equality")
})
