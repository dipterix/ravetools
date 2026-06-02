
# Helper: gain magnitude at normalized frequency f0 in [0,1]
fir_gain <- function(b, f0) {
  L <- length(b)
  Mod(sum(exp(-1i * 2 * pi * seq(0, L - 1) * (f0 / 2)) * b))
}

# ── Unity-gain tests (self-contained, no MATLAB needed) ──────────────────────

test_that("fir1: scale=TRUE produces unity gain at the correct reference frequency", {

  # Low-pass  -> DC gain = 1
  b <- fir1(20, 0.3, "low", hamming, scale = TRUE)$b
  expect_equal(fir_gain(b, 0), 1.0, tolerance = 1e-10)

  # High-pass -> Nyquist gain = 1
  b <- fir1(20, 0.3, "high", hamming, scale = TRUE)$b
  expect_equal(fir_gain(b, 1), 1.0, tolerance = 1e-10)

  # Band-pass -> gain at center of passband = 1
  # fir1 internally converts "pass" (2 cutoffs) to "DC-0" and uses mean(w)
  b <- fir1(20, c(0.2, 0.5), "pass", hamming, scale = TRUE)$b
  expect_equal(fir_gain(b, mean(c(0.2, 0.5))), 1.0, tolerance = 1e-10)

  # Band-stop -> DC gain = 1 (MATLAB fir1 convention: fband=TRUE for stop)
  b <- fir1(20, c(0.2, 0.5), "stop", hamming, scale = TRUE)$b
  expect_equal(fir_gain(b, 0), 1.0, tolerance = 1e-10)

  # scale=FALSE must not change the scaling -> gain differs from 1 noticeably
  b_unscaled <- fir1(20, 0.3, "low", hamming, scale = FALSE)$b
  b_scaled   <- fir1(20, 0.3, "low", hamming, scale = TRUE)$b
  expect_false(isTRUE(all.equal(b_unscaled, b_scaled)))
})


test_that("design_filter_fir: scale=TRUE produces unity gain at reference frequency (all methods × all types)", {

  fs   <- 1000
  nyq  <- fs / 2

  for (method in c("kaiser", "firls", "remez")) {

    # ── Low-pass: reference = DC (f0 = 0) ────────────────────────────────
    f <- design_filter_fir(
      sample_rate = fs, low_pass_freq = 100,
      filter_order = 60, scale = TRUE, method = method
    )
    expect_equal(
      fir_gain(f$b, 0), 1.0, tolerance = 1e-9,
      label = sprintf("low-pass %s DC gain", method)
    )

    # ── High-pass: reference = Nyquist (f0 = 1) ──────────────────────────
    f <- design_filter_fir(
      sample_rate = fs, high_pass_freq = 200,
      filter_order = 60, scale = TRUE, method = method
    )
    expect_equal(
      fir_gain(f$b, 1), 1.0, tolerance = 1e-9,
      label = sprintf("high-pass %s Nyquist gain", method)
    )

    # ── Band-pass: reference = mean(kaisprm$Wc) = mean of transition-band midpoints ────
    # high_pass_freq=100 <= low_pass_freq=300 -> ftype="pass"
    f <- design_filter_fir(
      sample_rate = fs, high_pass_freq = 100, low_pass_freq = 300,
      filter_order = 60, scale = TRUE, method = method
    )
    # mean(Wc) = (w_stop[1]+w_pass[1]+w_pass[2]+w_stop[2]) / 4
    f0_pass <- (f$parameters$stopband[[1]] + f$parameters$passband[[1]] +
                f$parameters$passband[[2]] + f$parameters$stopband[[2]]) / 4
    expect_equal(
      fir_gain(f$b, f0_pass), 1.0, tolerance = 1e-9,
      label = sprintf("band-pass %s transition-midpoint gain", method)
    )

    # ── Band-stop: reference = DC (f0 = 0) for both low- and high-freq notch ──
    # MATLAB fir1 always scales band-stop at DC (fband=TRUE -> b/sum(b))
    f_lo <- design_filter_fir(
      sample_rate = fs,
      low_pass_freq = 100, low_pass_trans_freq = 10,
      high_pass_freq = 200, high_pass_trans_freq = 10,
      filter_order = 60, scale = TRUE, method = method
    )
    expect_equal(
      fir_gain(f_lo$b, 0), 1.0, tolerance = 1e-9,
      label = sprintf("band-stop (low-freq notch) %s DC gain", method)
    )

    f_hi <- design_filter_fir(
      sample_rate = fs,
      low_pass_freq = 400, low_pass_trans_freq = 10,
      high_pass_freq = 450, high_pass_trans_freq = 10,
      filter_order = 60, scale = TRUE, method = method
    )
    expect_equal(
      fir_gain(f_hi$b, 0), 1.0, tolerance = 1e-9,
      label = sprintf("band-stop (high-freq notch) %s DC gain", method)
    )
  }
})


test_that("design_filter_fir: scale=FALSE does not normalise", {
  f_on  <- design_filter_fir(sample_rate = 1000, low_pass_freq = 100,
                              filter_order = 40, scale = TRUE)
  f_off <- design_filter_fir(sample_rate = 1000, low_pass_freq = 100,
                              filter_order = 40, scale = FALSE)
  # Scaled must have DC = 1; unscaled coefficients must differ
  expect_equal(fir_gain(f_on$b, 0), 1.0, tolerance = 1e-9)
  expect_false(isTRUE(all.equal(f_on$b, f_off$b)))
})


# ── MATLAB coefficient comparison ────────────────────────────────────────────
#
# Run the block below in MATLAB, paste the printed vectors back here, then
# remove the skip() calls.
#
# MATLAB code (paste into MATLAB command window or script):
# ──────────────────────────────────────────────────────────────────────────────
#   % fir1 defaults: Hamming window + scale=true, identical to ravetools defaults
#   b_low  = fir1(20, 0.3);
#   b_high = fir1(20, 0.3, 'high');
#   b_pass = fir1(20, [0.2 0.5], 'bandpass');
#   b_stop = fir1(20, [0.2 0.5], 'stop');
#
#   fprintf('b_low  <- c('); fprintf('%.15g, ', b_low);  fprintf(')\n');
#   fprintf('b_high <- c('); fprintf('%.15g, ', b_high); fprintf(')\n');
#   fprintf('b_pass <- c('); fprintf('%.15g, ', b_pass); fprintf(')\n');
#   fprintf('b_stop <- c('); fprintf('%.15g, ', b_stop); fprintf(')\n');
# ──────────────────────────────────────────────────────────────────────────────

test_that("fir1 low-pass matches MATLAB (hamming window, scale=TRUE)", {
  b_matlab <- c(
    0, 0.00292890317900091, 0.00634234732232293, 0.00378304197175637,
    -0.0123878480564484, -0.0343265728270622, -0.0318598610890541,
    0.0265312166066068, 0.137863164556997, 0.251347679146418,
    0.299555858378927, 0.251347679146418, 0.137863164556997,
    0.0265312166066068, -0.0318598610890541, -0.0343265728270622,
    -0.0123878480564484, 0.00378304197175637, 0.00634234732232293,
    0.00292890317900091, 0
  )

  b_r <- fir1(20, 0.3, "low", hamming, scale = TRUE)$b
  expect_equal(b_r, b_matlab, tolerance = 1e-7)
})


test_that("fir1 high-pass matches MATLAB (hamming window, scale=TRUE)", {
  b_matlab <- c(
    0, -0.00293014380119296, -0.00634503380813584, -0.00378464438929522,
    0.0123930952900616, 0.0343411128461149, 0.0318733562605525,
    -0.0265424546756356, -0.137921560513081, -0.251454144771811,
    0.699259735957564, -0.251454144771811, -0.137921560513081,
    -0.0265424546756356, 0.0318733562605525, 0.0343411128461149,
    0.0123930952900616, -0.00378464438929522, -0.00634503380813584,
    -0.00293014380119296, 0
  )

  b_r <- fir1(20, 0.3, "high", hamming, scale = TRUE)$b
  expect_equal(b_r, b_matlab, tolerance = 1e-7)
})


test_that("fir1 band-pass matches MATLAB (hamming window, scale=TRUE)", {
  b_matlab <- c(
    0, 0.0058661414341214, 0.00647237119566064, -0.00061145861614182,
    0.0126418101787922, 0.0350302987081212, -0.0325130171418959,
    -0.17094564692336, -0.140689484487979, 0.130693550172259,
    0.305697025206915, 0.130693550172259, -0.140689484487979,
    -0.17094564692336, -0.0325130171418959, 0.0350302987081212,
    0.0126418101787922, -0.00061145861614182, 0.00647237119566064,
    0.0058661414341214, 0
  )

  b_r <- fir1(20, c(0.2, 0.5), "pass", hamming, scale = TRUE)$b
  expect_equal(b_r, b_matlab, tolerance = 1e-7)
})


test_that("fir1 band-stop matches MATLAB (hamming window, scale=TRUE)", {
  b_matlab <- c(
    0, -0.00574321408319063, -0.00633674005647684, 0.000598645255139493,
    -0.0123768959666647, -0.0342962247225443, 0.0318316938030531,
    0.167363412201366, 0.137741279807586, -0.127954814312097,
    0.698345716147657, -0.127954814312097, 0.137741279807586,
    0.167363412201366, 0.0318316938030531, -0.0342962247225443,
    -0.0123768959666647, 0.000598645255139493, -0.00633674005647684,
    -0.00574321408319063, 0
  )

  b_r <- fir1(20, c(0.2, 0.5), "stop", hamming, scale = TRUE)$b
  expect_equal(b_r, b_matlab, tolerance = 1e-7)
})


# ── MATLAB comparison for design_filter_fir (kaiser method) ──────────────────
#
# The kaiser method calls fir1 with cutoffs = kaisprm$Wc and
# window = kaiser(n+1, beta).  Extract these from R, then verify in MATLAB.
#
# Run in R first:
#   f <- design_filter_fir(sample_rate=1000, low_pass_freq=100,
#                           filter_order=40, data_size=5000)
#   cat("n   =", f$parameters$order, "\n")
#   cat("Wc  =", f$parameters$cutoff, "\n")
#   cat("beta=", f$parameters$beta,  "\n")
#
# Then in MATLAB:
#   n = <n from R>; Wc = <Wc from R>; beta = <beta from R>;
#   b_lp_k = fir1(n, Wc, kaiser(n+1, beta), 'scale');
#   fprintf('b_lp_kaiser <- c('); fprintf('%.15g, ', b_lp_k); fprintf(')\n');

test_that("design_filter_fir kaiser low-pass matches MATLAB fir1", {
  # MATLAB run with: n=40; Wc=0.2013379; beta=3.395321;
  # b_lp_k = fir1(n, Wc, kaiser(n+1, beta), 'scale');
  # Note: Wc/beta were printed by R's cat() with ~7 sig figs; tolerance is loose accordingly.
  b_matlab <- c(
    0.000198963867774825, -0.00169241569059748, -0.00396307591421208,
    -0.00534352793206251, -0.0044412511715439, -0.000541195062692525,
    0.00567067137122843, 0.0118991132812478, 0.0148359902757939,
    0.0115083955960501, 0.000926301843132593, -0.0146875568050966,
    -0.0297746373755591, -0.0367515288539371, -0.0285015519433629,
    -0.00123060885531412, 0.0434159169471656, 0.0979264203726522,
    0.150565797367715, 0.188682526152012, 0.202594505059212,
    0.188682526152012, 0.150565797367715, 0.0979264203726522,
    0.0434159169471656, -0.00123060885531412, -0.0285015519433629,
    -0.0367515288539371, -0.0297746373755591, -0.0146875568050966,
    0.000926301843132593, 0.0115083955960501, 0.0148359902757939,
    0.0118991132812478, 0.00567067137122843, -0.000541195062692525,
    -0.0044412511715439, -0.00534352793206251, -0.00396307591421208,
    -0.00169241569059748, 0.000198963867774825
  )

  f <- design_filter_fir(
    sample_rate = 1000, low_pass_freq = 100,
    filter_order = 40, data_size = 5000
  )
  # Loose tolerance: Wc/beta were entered in MATLAB with only ~7 significant figures
  expect_equal(f$b, b_matlab, tolerance = 1e-4)
})

# ── Narrow band-stop (notch narrower than the auto transition heuristic) ─────

test_that("design_filter_fir: narrow band-stop with auto transitions succeeds", {
  # Notch is only 1 Hz wide (24-25 Hz) at sr=2000.  With data_size=NA the
  # auto suggested_trans_bandwidth is 0.01 (normalized) = 10 Hz, which is
  # 10x the notch.  The function must silently cap it rather than error.
  # Note: "remez" requires non-degenerate (non-zero-width) stop bands in its
  # internal grid and errors on near-zero transitions; that is a separate
  # gsignal constraint unrelated to this fix.
  for (method in c("kaiser", "firls")) {
    f <- design_filter_fir(
      sample_rate = 2000,
      high_pass_freq = 25, low_pass_freq = 24,
      filter_order = 200, method = method
    )
    expect_true(inherits(f, "ravetools-design_filter_fir"),
                label = sprintf("narrow notch returns filter (%s)", method))
    # Band-stop filters normalise at DC (f0 = 0)
    expect_equal(
      fir_gain(f$b, 0), 1.0, tolerance = 1e-6,
      label = sprintf("narrow notch DC gain = 1 (%s)", method)
    )
  }
})

test_that("design_filter_fir: narrow notch with one explicit transition succeeds", {
  # low_pass_trans_freq explicitly set to 0.3 Hz (< notch width 1 Hz);
  # high_pass_trans_freq is NA and will be capped to the remaining room.
  f <- design_filter_fir(
    sample_rate = 2000,
    high_pass_freq = 25, low_pass_freq = 24,
    low_pass_trans_freq = 0.3,
    filter_order = 200
  )
  expect_true(inherits(f, "ravetools-design_filter_fir"))
  expect_equal(fir_gain(f$b, 0), 1.0, tolerance = 1e-6)
})

test_that("design_filter_fir: user-supplied transitions that overlap throw an error", {
  # Combined explicit transition (1 + 1 = 2 Hz) exceeds 1 Hz notch -> error.
  expect_error(
    design_filter_fir(
      sample_rate = 2000,
      high_pass_freq = 25, low_pass_freq = 24,
      low_pass_trans_freq = 1, high_pass_trans_freq = 1,
      filter_order = 200
    ),
    regexp = "transition bandwidths overlap"
  )
  # A single explicit transition (1.5 Hz) that alone exceeds the notch -> error.
  expect_error(
    design_filter_fir(
      sample_rate = 2000,
      high_pass_freq = 25, low_pass_freq = 24,
      low_pass_trans_freq = 1.5,
      filter_order = 200
    ),
    regexp = "transition bandwidths overlap"
  )
})
