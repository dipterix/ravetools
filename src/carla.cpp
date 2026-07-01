// CARLA bootstrap inner loop (see R/carla.R).
//
// For one candidate subset size `ii`, `carla_zmin_boot` evaluates every
// bootstrap resampling entirely in C++: it trial-averages the resampled
// channels, re-references them by their own channel mean (the bootstrapped
// common average), correlates each unreferenced channel against every
// re-referenced channel over time, Fisher z-transforms, and returns - per
// bootstrap - the mean Fisher-z of the most globally anti-correlated channel.
//
// This replaces the R `M %*% W` resampling and the per-bootstrap `stats::cor`
// reduction, both of which dominated the runtime on large arrays. Bootstraps
// are independent and run in parallel via TinyParallel (honours
// `ravetools_threads()`); the random resampling indices are drawn in R and
// passed in as `ind` so the RNG stays reproducible and thread-safe.

#include <RcppEigen.h>
#include <cmath>
#include "utils.h"
#include "TinyParallel.h"

namespace {

struct CarlaZminWorker : public TinyParallel::Worker {
  const double* M;        // (ii*n_t) x n_tr, column-major (= `sub`)
  const int*    ind;      // n_resample x nboot, column-major, 1-based
  double*       out;      // nboot
  const int ii;
  const int n_t;
  const int n_tr;
  const int n_resample;

  CarlaZminWorker(const double* M, const int* ind, double* out,
                  int ii, int n_t, int n_tr, int n_resample)
    : M(M), ind(ind), out(out), ii(ii), n_t(n_t), n_tr(n_tr),
      n_resample(n_resample) {}

  void operator()(std::size_t begin, std::size_t end) {
    const Eigen::Index block = static_cast<Eigen::Index>(ii) *
      static_cast<Eigen::Index>(n_t);
    Eigen::Map<const Eigen::MatrixXd> Mmap(M, block, n_tr);
    const double inv = 1.0 / static_cast<double>(n_resample);
    const std::size_t chunk = end - begin;
    if (chunk == 0) return;

    // ---- 1. bootstrap weights for every column of `ind` in this chunk -----
    // Wc(j, c) = (number of times trial j is resampled in bootstrap begin+c)
    //            / n_resample. Then `M * Wc` trial-averages all bootstraps in
    // this chunk at once.
    Eigen::MatrixXd Wc =
      Eigen::MatrixXd::Zero(n_tr, static_cast<Eigen::Index>(chunk));
    for (std::size_t c = 0; c < chunk; c++) {
      const int* ind_b = ind + (begin + c) * static_cast<std::size_t>(n_resample);
      for (int k = 0; k < n_resample; k++) {
        const int j = ind_b[k];                  // 1-based trial index
        if (j >= 1 && j <= n_tr) {
          Wc(j - 1, static_cast<Eigen::Index>(c)) += inv;
        }
      }
    }

    // ---- 2. trial-averaged (resampled) means for the whole chunk ----------
    // One GEMM per thread reads the big `sub` matrix `M` only once (instead of
    // once per bootstrap), which is what lets this scale past memory bandwidth.
    // Column c is `Useg_m` (ii x n_t, column-major) for bootstrap begin+c.
    Eigen::MatrixXd Useg_chunk = Mmap * Wc;       // (ii*n_t) x chunk

    // Per-boot scratch, reused across the chunk to avoid reallocation.
    Eigen::MatrixXd    Ac(ii, n_t), B(ii, n_t), cross(ii, ii);
    Eigen::RowVectorXd cb(n_t);
    Eigen::VectorXd    am(ii), bm(ii), nrmA(ii), nrmB(ii);

    for (std::size_t c = 0; c < chunk; c++) {
      Eigen::Map<Eigen::MatrixXd> A(
          Useg_chunk.col(static_cast<Eigen::Index>(c)).data(), ii, n_t);

      // ---- 3. re-reference by channel mean (= bootstrapped CAR) -----------
      cb = A.colwise().mean();                    // length n_t
      B  = A;
      B.rowwise() -= cb;                          // Uref_m = Useg_m - cb

      // ---- 4. row-center over time + row norms ----------------------------
      Ac = A;
      am = Ac.rowwise().mean();
      Ac.colwise() -= am;
      bm = B.rowwise().mean();
      B.colwise()  -= bm;

      nrmA = Ac.rowwise().norm();                 // length ii
      nrmB = B.rowwise().norm();

      // ---- 5. cross products -> Pearson correlations over time ------------
      cross.noalias() = Ac * B.transpose();       // ii x ii ; sum_t Ac(j,t)B(l,t)

      // ---- 6. min over rows of the mean off-diagonal Fisher-z -------------
      // Mirrors R: z = atanh(cor); diag(z) <- NA; min over j of
      // rowMeans(z, na.rm = TRUE). Self (l == j), zero-variance rows / columns
      // and non-finite z are skipped (matches `na.rm = TRUE`).
      double best = R_PosInf;
      bool found = false;
      for (int j = 0; j < ii; j++) {
        if (!(nrmA[j] > 0.0)) continue;
        const double denomA = nrmA[j];
        double sum = 0.0;
        int cnt = 0;
        for (int l = 0; l < ii; l++) {
          if (l == j || !(nrmB[l] > 0.0)) continue;
          const double rr = cross(j, l) / (denomA * nrmB[l]);
          const double zz = std::atanh(rr);
          if (ISNAN(zz)) continue;
          sum += zz;
          cnt++;
        }
        if (cnt > 0) {
          const double m = sum / static_cast<double>(cnt);
          if (m < best) {
            best = m;
            found = true;
          }
        }
      }
      out[begin + c] = found ? best : NA_REAL;
    }
  }
};

} // anonymous namespace

// [[Rcpp::export(rng = false)]]
SEXP carla_zmin_boot(SEXP sub, SEXP ind) {

  // ---- dims of `sub` (ii x n_t x n_tr) -----------------------------------
  SEXP subDim = PROTECT(Rf_getAttrib(sub, R_DimSymbol));
  if (subDim == R_NilValue || Rf_length(subDim) != 3) {
    UNPROTECT(1);
    return make_error("C++ `carla_zmin_boot`: `sub` must be a 3-D array "
                      "(channels x time x trials).");
  }
  int ii, n_t, n_tr;
  if (TYPEOF(subDim) == REALSXP) {
    ii   = static_cast<int>(REAL(subDim)[0]);
    n_t  = static_cast<int>(REAL(subDim)[1]);
    n_tr = static_cast<int>(REAL(subDim)[2]);
  } else {
    ii   = INTEGER(subDim)[0];
    n_t  = INTEGER(subDim)[1];
    n_tr = INTEGER(subDim)[2];
  }
  UNPROTECT(1); // subDim

  // ---- coerce inputs to the expected storage types ------------------------
  SEXP sub_ = PROTECT(TYPEOF(sub) == REALSXP ? sub : Rf_coerceVector(sub, REALSXP));
  SEXP ind_ = PROTECT(TYPEOF(ind) == INTSXP  ? ind : Rf_coerceVector(ind, INTSXP));

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

  SEXP re = PROTECT(Rf_allocVector(REALSXP, nboot));

  CarlaZminWorker worker(REAL(sub_), INTEGER(ind_), REAL(re),
                         ii, n_t, n_tr, n_resample);
  TinyParallel::parallelFor(0, static_cast<std::size_t>(nboot), worker);

  UNPROTECT(3); // sub_, ind_, re
  return re;
}
