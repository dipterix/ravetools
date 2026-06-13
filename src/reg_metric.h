#ifndef RAVETOOLS_REG_METRIC_H
#define RAVETOOLS_REG_METRIC_H

// Similarity metrics + the per-level optimization problem for linear
// registration. Every metric reports a COST TO MINIMIZE and its gradient with
// respect to the transform's incremental parameters, so the optimizer always
// minimizes (information-theoretic / correlation metrics are negated).
//
// The expensive per-sample work (mapping each fixed point through the current
// transform and trilinearly sampling the moving image + its gradient) is done
// once per metric evaluation, in parallel, into per-sample caches; the metric
// math (histogram, correlation sums, gradient accumulation) then reads those
// caches serially.

#include <RcppEigen.h>
#include <vector>
#include <random>
#include <cmath>
#include "reg_core.h"
#include "reg_transform.h"
#include "TinyParallel.h"

namespace ravereg {

using Eigen::VectorXd;

enum class Metric { MEANSQUARES, MATTES, CC };

// Parallel worker: map each cached fixed RAS point through the transform and
// sample the moving image (value, and optionally RAS gradient). When a moving
// mask is supplied, a sample whose mapped point lands outside it (mask <= 0.5)
// is marked invalid so it drops out of the metric.
struct MovingSampler : public TinyParallel::Worker {
  const RegImage3* mov;
  const RegImage3* movMask = nullptr;   // optional moving mask (nullptr = none)
  Matrix3d A; Vector3d t; Vector3d center;
  const double *px, *py, *pz;
  double *mw, *gx, *gy, *gz;
  unsigned char* valid;
  bool needGrad;

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t s = begin; s < end; ++s) {
      Vector3d p(px[s], py[s], pz[s]);
      Vector3d q = A * (p - center) + center + t;
      double val;
      Vector3d g;
      bool ok = mov->sampleRas(q, val, needGrad ? &g : nullptr);
      if (ok && movMask) {
        double mv; movMask->sampleRas(q, mv);
        if (mv <= 0.5) ok = false;        // mapped outside the moving mask
      }
      valid[s] = ok ? 1 : 0;
      mw[s] = val;
      if (needGrad) { gx[s] = g[0]; gy[s] = g[1]; gz[s] = g[2]; }
    }
  }
};

struct LinearProblem {
  const FixedLevel* fixed = nullptr;
  const RegImage3* moving = nullptr;
  LinearTransform xform;
  Metric metric = Metric::MEANSQUARES;

  // Optional masks (nullptr = none). fixedMaskLevel is parallel to fixed->data
  // (same level grid): only its non-zero voxels are sampled. movingMask is a
  // full-moving-grid image: a sample mapping outside it is dropped.
  const reg_real* fixedMaskLevel = nullptr;
  const RegImage3* movingMask = nullptr;

  // sampled fixed-voxel linear indices for this level
  std::vector<std::size_t> sampleIdx;

  // per-sample caches (sized to sampleIdx)
  std::vector<double> sPx, sPy, sPz;     // fixed RAS point (constant per level)
  std::vector<double> sF;                // fixed intensity
  std::vector<double> sMw, sGx, sGy, sGz; // warped moving value + RAS gradient
  std::vector<unsigned char> sValid;

  // Mattes MI configuration
  int numberOfBins = 32;
  int ccRadius = 2;   // (linear CC is global; field reserved for parity)

  double fixedMin = 0.0, fixedMax = 1.0;
  double movingMin = 0.0, movingMax = 1.0;

  // ----- helpers ---------------------------------------------------------
  inline void decode(std::size_t lin, double& i, double& j, double& k) const {
    const std::size_t syf = static_cast<std::size_t>(fixed->nx);
    const std::size_t szf = syf * static_cast<std::size_t>(fixed->ny);
    const std::size_t kk = lin / szf;
    const std::size_t r = lin - kk * szf;
    const std::size_t jj = r / syf;
    const std::size_t ii = r - jj * syf;
    i = static_cast<double>(ii);
    j = static_cast<double>(jj);
    k = static_cast<double>(kk);
  }

  inline Vector3d fixedRas(std::size_t lin) const {
    double i, j, k;
    decode(lin, i, j, k);
    Vector4d vh(i, j, k, 1.0);
    return (fixed->vox2ras * vh).head<3>();
  }

  void buildSamples(double rate, unsigned int seed) {
    sampleIdx.clear();
    const std::size_t n = fixed->size();
    // A fixed mask restricts sampling to its non-zero voxels; the sampling rate
    // then thins within that set (so masked registration is also faster).
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    sampleIdx.reserve(static_cast<std::size_t>(n * (rate >= 1.0 ? 1.0 : rate) * 1.1));
    for (std::size_t s = 0; s < n; ++s) {
      if (fixedMaskLevel && fixedMaskLevel[s] <= reg_real(0)) continue;
      if (rate >= 1.0 || unif(rng) < rate) sampleIdx.push_back(s);
    }
    if (sampleIdx.empty()) sampleIdx.push_back(0);
    // cache constant per-sample fixed RAS points + fixed intensities
    const std::size_t m = sampleIdx.size();
    sPx.resize(m); sPy.resize(m); sPz.resize(m); sF.resize(m);
    for (std::size_t s = 0; s < m; ++s) {
      const Vector3d p = fixedRas(sampleIdx[s]);
      sPx[s] = p[0]; sPy[s] = p[1]; sPz[s] = p[2];
      sF[s] = fixed->data[sampleIdx[s]];
    }
    sMw.resize(m); sGx.resize(m); sGy.resize(m); sGz.resize(m); sValid.resize(m);
  }

  void prepareRanges() {
    fixedMin = std::numeric_limits<double>::max();
    fixedMax = -std::numeric_limits<double>::max();
    const std::size_t n = fixed->size();
    for (std::size_t s = 0; s < n; ++s) {
      const double v = fixed->data[s];
      if (v < fixedMin) fixedMin = v;
      if (v > fixedMax) fixedMax = v;
    }
    movingMin = std::numeric_limits<double>::max();
    movingMax = -std::numeric_limits<double>::max();
    const std::size_t m = moving->size();
    for (std::size_t s = 0; s < m; ++s) {
      const double v = moving->data[s];
      if (v < movingMin) movingMin = v;
      if (v > movingMax) movingMax = v;
    }
    if (fixedMax <= fixedMin) fixedMax = fixedMin + 1.0;
    if (movingMax <= movingMin) movingMax = movingMin + 1.0;
  }

  // Parallel: refresh the warped-moving caches for the current transform.
  void sampleMoving(bool needGrad) {
    MovingSampler w;
    w.mov = moving;
    w.movMask = movingMask;
    w.A = xform.A; w.t = xform.t; w.center = xform.center;
    w.px = sPx.data(); w.py = sPy.data(); w.pz = sPz.data();
    w.mw = sMw.data(); w.gx = sGx.data(); w.gy = sGy.data(); w.gz = sGz.data();
    w.valid = sValid.data(); w.needGrad = needGrad;
    TinyParallel::parallelFor(0, sampleIdx.size(), w, 256);
  }

  // ----- parameter scales (physical-shift) -------------------------------
  void estimateScales(VectorXd& scales) const {
    const int np = xform.nParams();
    VectorXd acc = VectorXd::Zero(np);
    long count = 0;
    for (std::size_t s = 0; s < sampleIdx.size(); ++s) {
      const Vector3d p(sPx[s], sPy[s], sPz[s]);
      const Vector3d yv = p - xform.center;
      if (xform.mode == LinearMode::RIGID) {
        const Vector3d w = xform.A * yv;
        const double wn2 = w.squaredNorm();
        acc[0] += wn2 - w[0] * w[0];
        acc[1] += wn2 - w[1] * w[1];
        acc[2] += wn2 - w[2] * w[2];
        acc[3] += 1.0; acc[4] += 1.0; acc[5] += 1.0;
      } else {
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            acc[a * 3 + b] += yv[b] * yv[b];
        acc[9] += 1.0; acc[10] += 1.0; acc[11] += 1.0;
      }
      ++count;
    }
    scales = VectorXd::Ones(np);
    if (count > 0) {
      acc /= static_cast<double>(count);
      for (int j = 0; j < np; ++j)
        scales[j] = (acc[j] > 1e-12) ? acc[j] : 1.0;
    }
  }

  // ----- metric dispatch -------------------------------------------------
  double evaluate(VectorXd* grad) {
    sampleMoving(grad != nullptr);
    switch (metric) {
      case Metric::MEANSQUARES: return evalMeanSquares(grad);
      case Metric::MATTES:      return evalMattes(grad);
      case Metric::CC:          return evalCC(grad);
    }
    return 0.0;
  }

  // Accumulate d(cost)/d(params) for a single sample.
  inline void accumulateLinearGrad(VectorXd& grad, double rf,
                                   const Vector3d& yv, const Vector3d& gM) const {
    if (xform.mode == LinearMode::RIGID) {
      const Vector3d w = xform.A * yv;
      const Vector3d rot = w.cross(gM);     // (-[w]_x)^T gM = w x gM
      grad[0] += rf * rot[0]; grad[1] += rf * rot[1]; grad[2] += rf * rot[2];
      grad[3] += rf * gM[0];  grad[4] += rf * gM[1];  grad[5] += rf * gM[2];
    } else {
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          grad[a * 3 + b] += rf * gM[a] * yv[b];
      grad[9]  += rf * gM[0];
      grad[10] += rf * gM[1];
      grad[11] += rf * gM[2];
    }
  }

  // ----- mean squares ----------------------------------------------------
  double evalMeanSquares(VectorXd* grad) {
    const int np = xform.nParams();
    if (grad) grad->setZero(np);
    double cost = 0.0;
    long count = 0;
    for (std::size_t s = 0; s < sampleIdx.size(); ++s) {
      if (!sValid[s]) continue;
      const double diff = sMw[s] - sF[s];
      cost += diff * diff;
      ++count;
      if (grad) {
        const Vector3d yv(sPx[s] - xform.center[0],
                          sPy[s] - xform.center[1],
                          sPz[s] - xform.center[2]);
        const Vector3d gM(sGx[s], sGy[s], sGz[s]);
        accumulateLinearGrad(*grad, 2.0 * diff, yv, gM);
      }
    }
    if (count > 0) {
      cost /= static_cast<double>(count);
      if (grad) *grad /= static_cast<double>(count);
    }
    return cost;
  }

  double evalMattes(VectorXd* grad);
  double evalCC(VectorXd* grad);
};

} // namespace ravereg

#endif // RAVETOOLS_REG_METRIC_H
