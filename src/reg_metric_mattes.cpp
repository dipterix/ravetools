// Mattes mutual information metric for linear registration.
//
// Follows Mattes et al. 2003 (the ITK MattesMutualInformationImageToImageMetric
// formulation): the fixed image is binned with a zero-order B-spline (one bin
// per sample) and the moving image with a cubic B-spline Parzen window, giving a
// differentiable joint histogram J[f][m].
//
//   MI   = sum_{f,m} J log( J / (Jf * Jm) )          (Jf is constant in mu)
//   cost = -MI
//   d(MI)/dmu = -(1/(N*dm)) sum_samples [ sum_m b3'(m - m*) log(J/Jm) ] dM/dmu
//
// so the per-sample residual factor handed to the shared linear-gradient
// accumulator is  bracket / (N * movingBinSize).
//
// Performance notes:
//   Phase 1 (histogram): parallelized via parallelReduce with per-thread
//     accumulators. fbin/mstar/mbase are cached for Phase 2 reuse.
//   Phase 2 (gradient): parallelized via parallelReduce with per-thread
//     gradient accumulators. log() calls replaced by a precomputed B² table.

#include "reg_metric.h"

namespace ravereg {

// Cubic B-spline kernel and its derivative.
static inline double bspline3(double t) {
  const double a = std::fabs(t);
  if (a < 1.0) return (4.0 - 6.0 * a * a + 3.0 * a * a * a) / 6.0;
  if (a < 2.0) { const double b = 2.0 - a; return b * b * b / 6.0; }
  return 0.0;
}

static inline double bspline3deriv(double t) {
  const double a = std::fabs(t);
  if (a < 1.0) return -2.0 * t + 1.5 * t * a;           // d/dt (2/3 - t^2 + |t|^3/2)
  if (a < 2.0) {
    const double s = (t >= 0.0) ? 1.0 : -1.0;
    const double b = 2.0 - a;
    return -0.5 * s * b * b;
  }
  return 0.0;
}

// ---- Phase 1: parallel histogram build + bin-state cache fill ---------------
//
// Each thread accumulates into its own J histogram to avoid write contention.
// fbin_c / mstar_c / mbase_c are written at sample-disjoint indices (safe).
// join() reduces per-thread histograms by element-wise addition (O(B²)).

struct MattesHistReducer : public TinyParallel::Worker {
  const double*        sF;
  const double*        sMw;
  const unsigned char* sValid;
  int    bins;
  double fixedMin, fixedBinSize;
  double movingMin, movingBinSize;
  double pad_d;

  // shared output caches — each sample index is unique to one thread's range
  int*    fbin_c;
  double* mstar_c;
  int*    mbase_c;

  std::vector<double> J;   // per-thread accumulator (bins × bins)

  MattesHistReducer(const double* sF_, const double* sMw_,
                    const unsigned char* sValid_,
                    int bins_,
                    double fixedMin_, double fixedBinSize_,
                    double movingMin_, double movingBinSize_, int pad_,
                    int* fbin_c_, double* mstar_c_, int* mbase_c_)
    : sF(sF_), sMw(sMw_), sValid(sValid_)
    , bins(bins_)
    , fixedMin(fixedMin_), fixedBinSize(fixedBinSize_)
    , movingMin(movingMin_), movingBinSize(movingBinSize_)
    , pad_d(static_cast<double>(pad_))
    , fbin_c(fbin_c_), mstar_c(mstar_c_), mbase_c(mbase_c_)
    , J(static_cast<std::size_t>(bins_) * bins_, 0.0)
  {}

  MattesHistReducer(MattesHistReducer& o, TinyParallel::Split)
    : sF(o.sF), sMw(o.sMw), sValid(o.sValid)
    , bins(o.bins)
    , fixedMin(o.fixedMin), fixedBinSize(o.fixedBinSize)
    , movingMin(o.movingMin), movingBinSize(o.movingBinSize)
    , pad_d(o.pad_d)
    , fbin_c(o.fbin_c), mstar_c(o.mstar_c), mbase_c(o.mbase_c)
    , J(static_cast<std::size_t>(o.bins) * o.bins, 0.0)
  {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t s = begin; s < end; ++s) {
      if (!sValid[s]) { fbin_c[s] = -1; continue; }

      int fb = static_cast<int>(std::floor((sF[s] - fixedMin) / fixedBinSize))
               + static_cast<int>(pad_d);
      if (fb < 0) fb = 0; else if (fb > bins - 1) fb = bins - 1;

      const double ms = (sMw[s] - movingMin) / movingBinSize + pad_d;
      const int mb    = static_cast<int>(std::floor(ms));

      fbin_c[s]  = fb;
      mstar_c[s] = ms;
      mbase_c[s] = mb;

      for (int k = mb - 1; k <= mb + 2; ++k) {
        if (k < 0 || k > bins - 1) continue;
        const double w = bspline3(static_cast<double>(k) - ms);
        if (w != 0.0) J[static_cast<std::size_t>(fb) * bins + k] += w;
      }
    }
  }

  void join(const MattesHistReducer& rhs) {
    const std::size_t n = J.size();
    for (std::size_t i = 0; i < n; ++i) J[i] += rhs.J[i];
  }
};

// ---- Phase 2: parallel gradient accumulation --------------------------------
//
// Uses cached fbin/mstar/mbase (no recomputation) and a precomputed B² logRatio
// table (replaces O(4N) log() calls with O(B²) log() calls + O(4N) lookups).
// Per-thread VectorXd accumulators are reduced by vector addition (O(np)).

struct MattesGradReducer : public TinyParallel::Worker {
  const int*    fbin_c;
  const double* mstar_c;
  const int*    mbase_c;
  const double* logRatio;          // B × B table: log(J[f,m] / Jm[m])
  const double* sPx, *sPy, *sPz;
  const double* sGx, *sGy, *sGz;
  const LinearTransform* xform;
  double scale;
  int bins, np;

  VectorXd localGrad;

  MattesGradReducer(const int* fbin_c_, const double* mstar_c_, const int* mbase_c_,
                    const double* logRatio_,
                    const double* sPx_, const double* sPy_, const double* sPz_,
                    const double* sGx_, const double* sGy_, const double* sGz_,
                    const LinearTransform* xform_, double scale_, int bins_, int np_)
    : fbin_c(fbin_c_), mstar_c(mstar_c_), mbase_c(mbase_c_)
    , logRatio(logRatio_)
    , sPx(sPx_), sPy(sPy_), sPz(sPz_)
    , sGx(sGx_), sGy(sGy_), sGz(sGz_)
    , xform(xform_), scale(scale_), bins(bins_), np(np_)
    , localGrad(VectorXd::Zero(np_))
  {}

  MattesGradReducer(MattesGradReducer& o, TinyParallel::Split)
    : fbin_c(o.fbin_c), mstar_c(o.mstar_c), mbase_c(o.mbase_c)
    , logRatio(o.logRatio)
    , sPx(o.sPx), sPy(o.sPy), sPz(o.sPz)
    , sGx(o.sGx), sGy(o.sGy), sGz(o.sGz)
    , xform(o.xform), scale(o.scale), bins(o.bins), np(o.np)
    , localGrad(VectorXd::Zero(o.np))
  {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t s = begin; s < end; ++s) {
      if (fbin_c[s] < 0) continue;
      const int    fb = fbin_c[s];
      const double ms = mstar_c[s];
      const int    mb = mbase_c[s];

      double bracket = 0.0;
      for (int k = mb - 1; k <= mb + 2; ++k) {
        if (k < 0 || k > bins - 1) continue;
        bracket += bspline3deriv(static_cast<double>(k) - ms)
                   * logRatio[static_cast<std::size_t>(fb) * bins + k];
      }

      const double rf = scale * bracket;
      const Vector3d yv(sPx[s] - xform->center[0],
                        sPy[s] - xform->center[1],
                        sPz[s] - xform->center[2]);
      const Vector3d gM(sGx[s], sGy[s], sGz[s]);

      // inline accumulateLinearGrad to write into thread-local localGrad
      if (xform->mode == LinearMode::RIGID) {
        const Vector3d w   = xform->A * yv;
        const Vector3d rot = w.cross(gM);
        localGrad[0] += rf * rot[0]; localGrad[1] += rf * rot[1]; localGrad[2] += rf * rot[2];
        localGrad[3] += rf * gM[0];  localGrad[4] += rf * gM[1];  localGrad[5] += rf * gM[2];
      } else {
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            localGrad[a * 3 + b] += rf * gM[a] * yv[b];
        localGrad[9]  += rf * gM[0];
        localGrad[10] += rf * gM[1];
        localGrad[11] += rf * gM[2];
      }
    }
  }

  void join(const MattesGradReducer& rhs) {
    localGrad += rhs.localGrad;
  }
};

// ---- evalMattes -------------------------------------------------------------

double LinearProblem::evalMattes(VectorXd* grad) {
  const int np = xform.nParams();
  if (grad) grad->setZero(np);

  const int bins = numberOfBins;
  const int pad  = 2;                     // cubic B-spline padding
  const double eps = 1e-16;

  const double fRange = fixedMax - fixedMin;
  const double mRange = movingMax - movingMin;
  const double fixedBinSize  = fRange / static_cast<double>(bins - 2 * pad);
  const double movingBinSize = mRange / static_cast<double>(bins - 2 * pad);
  if (fixedBinSize <= 0.0 || movingBinSize <= 0.0) return 0.0;

  const std::size_t ns = sampleIdx.size();

  // Per-sample bin-state caches (computed once in Phase 1, reused in Phase 2)
  std::vector<int>    fbin_c(ns);
  std::vector<double> mstar_c(ns);
  std::vector<int>    mbase_c(ns);

  // ---- Phase 1 (parallel): histogram accumulation + cache fill ------------
  MattesHistReducer histR(
    sF.data(), sMw.data(), sValid.data(),
    bins, fixedMin, fixedBinSize, movingMin, movingBinSize, pad,
    fbin_c.data(), mstar_c.data(), mbase_c.data());
  TinyParallel::parallelReduce(0, ns, histR, 256);

  std::vector<double> J = std::move(histR.J);

  // Count valid samples from cache (fbin_c[s] == -1 means invalid)
  long nvalid = 0;
  for (std::size_t s = 0; s < ns; ++s)
    if (fbin_c[s] >= 0) ++nvalid;

  const double N = static_cast<double>(nvalid);
  if (N < 1.0) return 0.0;

  // Normalize and build marginals
  for (double& v : J) v /= N;
  std::vector<double> Jm(bins, 0.0);
  std::vector<double> Jf(bins, 0.0);
  for (int f = 0; f < bins; ++f)
    for (int m = 0; m < bins; ++m) {
      const double v = J[static_cast<std::size_t>(f) * bins + m];
      Jf[f] += v;
      Jm[m] += v;
    }

  // ---- Mutual information (unchanged from serial version) -----------------
  double mi = 0.0;
  for (int f = 0; f < bins; ++f) {
    if (Jf[f] < eps) continue;
    for (int m = 0; m < bins; ++m) {
      const double v = J[static_cast<std::size_t>(f) * bins + m];
      if (v < eps || Jm[m] < eps) continue;
      mi += v * std::log(v / (Jf[f] * Jm[m]));
    }
  }

  // ---- Phase 2 (parallel): gradient accumulation --------------------------
  if (grad) {
    // Precompute logRatio[f*bins+m] = log(J[f,m] / Jm[m])  (B² log calls)
    std::vector<double> logRatio(static_cast<std::size_t>(bins) * bins, 0.0);
    for (int f = 0; f < bins; ++f)
      for (int m = 0; m < bins; ++m) {
        const double Jval = J[static_cast<std::size_t>(f) * bins + m];
        if (Jval > eps && Jm[m] > eps)
          logRatio[static_cast<std::size_t>(f) * bins + m] = std::log(Jval / Jm[m]);
      }

    const double scale = 1.0 / (N * movingBinSize);
    MattesGradReducer gradR(
      fbin_c.data(), mstar_c.data(), mbase_c.data(),
      logRatio.data(),
      sPx.data(), sPy.data(), sPz.data(),
      sGx.data(), sGy.data(), sGz.data(),
      &xform, scale, bins, np);
    TinyParallel::parallelReduce(0, ns, gradR, 256);
    *grad = gradR.localGrad;
  }

  return -mi;
}

} // namespace ravereg
