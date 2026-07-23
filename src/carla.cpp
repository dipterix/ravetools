// CARLA bootstrap evaluation (see R/carla.R).
//
// `carla_zmin_all` evaluates the CARLA anti-correlation diagnostic for every
// candidate subset size ii = 2..n_good against ONE shared set of bootstrap
// trial-resamples. The work is split into two phases:
//
//   Phase A (parallel over channels): trial-average every channel for every
//     bootstrap at once via a single GEMM per channel, `Sc * W`, where `W` is
//     the n_tr x nboot bootstrap-weight matrix. This reads the (large) input
//     `x_ord` exactly once, instead of once per bootstrap, so the bandwidth-
//     bound averaging is paid a single time.
//
//   Phase B (parallel over bootstraps): for each bootstrap, center the channel
//     averages over time, form the channel x channel Gram matrix G = Ac^T Ac
//     once, and read every subset size out of G in O(ii^2) scalar work -- so
//     the cross-products are paid once per bootstrap rather than once per
//     subset size (the old per-subset-size loop re-trial-averaged everything).
//
// Bootstraps / channels are independent and run in parallel via TinyParallel
// (honours `ravetools_threads()`); the resampling indices are drawn in R and
// passed in as `ind` so the RNG stays reproducible and thread-safe.
//
// Math. With Ac_c = time-centered trial-average of channel c, G = Ac^T Ac,
// and running sums p[j] = sum_{c<=ii} G[j,c], Q = sum_{j,l<=ii} G[j,l]:
//   cross(j,l) = G[j,l] - p[j]/ii            (= <Ac_j, Ac_l - c0>, c0 = mean CAR)
//   ||B_l||^2  = G[l,l] - 2 p[l]/ii + Q/ii^2 (= ||Ac_l - c0||^2, a true norm)
//   r(j,l) = cross / (||Ac_j|| * ||B_l||); z = atanh(r); drop l==j and any
//   zero-variance / non-finite entry (matches R `na.rm = TRUE`);
//   zmin[ii] = min_j mean_{l != j} z(j,l).
// p and Q update incrementally as ii grows, so each subset size costs O(ii^2).

#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include "utils.h"
#include "TinyParallel.h"

namespace {

// ---- Phase A: batched trial-averaging, one GEMM per channel ---------------
struct TrialAvgWorker : public TinyParallel::Worker {
  const double* M;        // n_t x n_tr x n_good, column-major (= `x_ord`)
  const double* W;        // n_tr x nboot, column-major (bootstrap weights)
  double*       U;        // n_t x nboot x n_good, column-major (output)
  const int n_t;
  const int n_tr;
  const int nboot;

  TrialAvgWorker(const double* M, const double* W, double* U,
                 int n_t, int n_tr, int nboot)
    : M(M), W(W), U(U), n_t(n_t), n_tr(n_tr), nboot(nboot) {}

  void operator()(std::size_t begin, std::size_t end) {
    Eigen::Map<const Eigen::MatrixXd> Wm(W, n_tr, nboot);
    for (std::size_t c = begin; c < end; c++) {
      Eigen::Map<const Eigen::MatrixXd> Sc(
          M + static_cast<std::size_t>(n_t) * n_tr * c, n_t, n_tr);
      Eigen::Map<Eigen::MatrixXd> Uc(
          U + static_cast<std::size_t>(n_t) * nboot * c, n_t, nboot);
      Uc.noalias() = Sc * Wm;                     // n_t x nboot trial-averages
    }
  }
};

// ---- Phase B: Gram matrix + subset-size sweep, per bootstrap --------------
struct CarlaZminWorker : public TinyParallel::Worker {
  const double* U;        // n_t x nboot x n_good, column-major (from Phase A)
  double*       out;      // n_good x nboot, column-major
  const int n_good;
  const int n_t;
  const int nboot;

  CarlaZminWorker(const double* U, double* out, int n_good, int n_t, int nboot)
    : U(U), out(out), n_good(n_good), n_t(n_t), nboot(nboot) {}

  void operator()(std::size_t begin, std::size_t end) {
    // Per-thread scratch, reused across every bootstrap in this range.
    Eigen::MatrixXd Ac(n_t, n_good);
    Eigen::MatrixXd G(n_good, n_good);
    Eigen::VectorXd diagG(n_good), normA(n_good), normB(n_good), p(n_good);

    for (std::size_t b = begin; b < end; b++) {

      // ---- 1. gather this bootstrap's channel averages, center over time --
      for (int c = 0; c < n_good; c++) {
        Eigen::Map<const Eigen::VectorXd> uc(
            U + static_cast<std::size_t>(n_t) * b +
                static_cast<std::size_t>(n_t) * nboot * c, n_t);
        Ac.col(c) = uc.array() - uc.mean();
      }

      // ---- 2. Gram matrix of the centered channel averages ----------------
      G.noalias() = Ac.transpose() * Ac;          // n_good x n_good, symmetric
      for (int c = 0; c < n_good; c++) {
        diagG[c] = G(c, c);
        normA[c] = diagG[c] > 0.0 ? std::sqrt(diagG[c]) : 0.0;
      }

      // ---- 3. sweep subset sizes ii = 2..n_good straight from the Gram ----
      double* out_b = out + static_cast<std::size_t>(n_good) * b;
      out_b[0] = NA_REAL;                          // subset size 1 is unused

      // Running p / Q for the currently-included channel prefix. Seed with
      // channel 0 (subset size 1, not evaluated).
      p.setZero();
      double Q = G(0, 0);
      for (int j = 0; j < n_good; j++) p[j] += G(j, 0);

      for (int m = 2; m <= n_good; m++) {
        const int mi = m - 1;                      // 0-based channel just added
        // Extend the prefix to include channel mi (uses p[mi] pre-update).
        Q += 2.0 * p[mi] + G(mi, mi);
        for (int j = 0; j < n_good; j++) p[j] += G(j, mi);

        const double invm  = 1.0 / static_cast<double>(m);
        const double invm2 = invm * invm;
        for (int l = 0; l < m; l++) {
          const double nb2 = diagG[l] - 2.0 * p[l] * invm + Q * invm2;
          normB[l] = nb2 > 0.0 ? std::sqrt(nb2) : 0.0;
        }

        double best = R_PosInf;
        bool found = false;
        for (int j = 0; j < m; j++) {
          if (!(normA[j] > 0.0)) continue;
          const double denomA = normA[j];
          const double pj = p[j] * invm;
          double sum = 0.0;
          int cnt = 0;
          for (int l = 0; l < m; l++) {
            if (l == j || !(normB[l] > 0.0)) continue;
            double rr = (G(j, l) - pj) / (denomA * normB[l]);
            // Clamp for FP safety: Gram cancellation can push |r| just past 1.
            if (rr >  1.0 - 1e-12) rr =  1.0 - 1e-12;
            if (rr < -1.0 + 1e-12) rr = -1.0 + 1e-12;
            const double zz = std::atanh(rr);
            if (ISNAN(zz)) continue;
            sum += zz;
            cnt++;
          }
          if (cnt > 0) {
            const double mmean = sum / static_cast<double>(cnt);
            if (mmean < best) {
              best = mmean;
              found = true;
            }
          }
        }
        out_b[m - 1] = found ? best : NA_REAL;
      }
    }
  }
};

} // anonymous namespace

// [[Rcpp::export(rng = false)]]
SEXP carla_zmin_all(SEXP x_ord, SEXP ind) {

  // ---- dims of `x_ord` (n_t x n_tr x n_good) -----------------------------
  SEXP subDim = PROTECT(Rf_getAttrib(x_ord, R_DimSymbol));
  if (subDim == R_NilValue || Rf_length(subDim) != 3) {
    UNPROTECT(1);
    return make_error("C++ `carla_zmin_all`: `x_ord` must be a 3-D array "
                      "(time x trials x channels).");
  }
  int n_good, n_t, n_tr;
  if (TYPEOF(subDim) == REALSXP) {
    n_t    = static_cast<int>(REAL(subDim)[0]);
    n_tr   = static_cast<int>(REAL(subDim)[1]);
    n_good = static_cast<int>(REAL(subDim)[2]);
  } else {
    n_t    = INTEGER(subDim)[0];
    n_tr   = INTEGER(subDim)[1];
    n_good = INTEGER(subDim)[2];
  }
  UNPROTECT(1); // subDim

  // ---- coerce inputs to the expected storage types ------------------------
  SEXP x_   = PROTECT(TYPEOF(x_ord) == REALSXP ? x_ord : Rf_coerceVector(x_ord, REALSXP));
  SEXP ind_ = PROTECT(TYPEOF(ind)   == INTSXP  ? ind   : Rf_coerceVector(ind, INTSXP));

  // ---- dims of `ind` (n_resample x nboot) --------------------------------
  int n_resample, nboot;
  SEXP indDim = PROTECT(Rf_getAttrib(ind_, R_DimSymbol));
  if (indDim != R_NilValue && Rf_length(indDim) == 2) {
    n_resample = INTEGER(indDim)[0];
    nboot      = INTEGER(indDim)[1];
  } else {
    n_resample = static_cast<int>(Rf_xlength(ind_));
    nboot      = 1;
  }
  UNPROTECT(1); // indDim

  // ---- bootstrap-weight matrix W (n_tr x nboot) --------------------------
  // W(j, b) = (number of times trial j is resampled in bootstrap b) / n_resample.
  const int* ind_p = INTEGER(ind_);
  const double inv = 1.0 / static_cast<double>(n_resample);
  std::vector<double> W(static_cast<std::size_t>(n_tr) * nboot, 0.0);
  for (int b = 0; b < nboot; b++) {
    const int* ind_b = ind_p + static_cast<std::size_t>(b) * n_resample;
    double* Wb = W.data() + static_cast<std::size_t>(b) * n_tr;
    for (int k = 0; k < n_resample; k++) {
      const int j = ind_b[k];                      // 1-based trial index
      if (j >= 1 && j <= n_tr) Wb[j - 1] += inv;
    }
  }

  // ---- Phase A: batched trial-averages, U (n_t x nboot x n_good) ----------
  std::vector<double> U(static_cast<std::size_t>(n_t) * nboot * n_good);
  TrialAvgWorker avg(REAL(x_), W.data(), U.data(), n_t, n_tr, nboot);
  TinyParallel::parallelFor(0, static_cast<std::size_t>(n_good), avg);

  // ---- Phase B: Gram + subset sweep, out (n_good x nboot) ----------------
  SEXP re = PROTECT(Rf_allocMatrix(REALSXP, n_good, nboot));
  CarlaZminWorker worker(U.data(), REAL(re), n_good, n_t, nboot);
  TinyParallel::parallelFor(0, static_cast<std::size_t>(nboot), worker);

  UNPROTECT(3); // x_, ind_, re
  return re;
}
