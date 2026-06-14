// Native linear (rigid / affine) 3D registration entry point.
//
// Multiresolution, physical-shift-scaled regular-step gradient descent over a
// similarity metric (mean squares, Mattes MI, or global normalized CC), working
// entirely in RAS physical space. Each image carries its own vox2ras.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include <iomanip>
#include "reg_core.h"
#include "reg_transform.h"
#include "reg_metric.h"

using namespace ravereg;

/*** R
# DIPSAUS DEBUG START
# devtools::load_all()
# x <- array(0, c(40, 40, 40)); x[12:28, 14:26, 10:30] <- 1
# v2r <- diag(c(1, 1, 1, 1)); v2r[1:3, 4] <- -20
# # shift the moving image by a known translation in RAS
# mov2ras <- v2r; mov2ras[1, 4] <- mov2ras[1, 4] - 3
# res <- register_linear_cpp(
#   as.double(x), dim(x), v2r, as.double(x), dim(x), mov2ras,
#   "rigid", "meansquares", c(4L, 2L, 1L), c(2, 1, 0), c(200L, 200L, 100L),
#   0.5, 0, 32L, 1L, NULL)
# res$transform   # should recover ~ +3 mm translation on x
*/

namespace {

Matrix4d asMatrix4(const Rcpp::NumericMatrix& m) {
  Matrix4d M = Matrix4d::Identity();
  const int nr = std::min(4, (int)m.nrow());
  const int nc = std::min(4, (int)m.ncol());
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j < nc; ++j)
      M(i, j) = m(i, j);
  return M;
}

Rcpp::NumericMatrix asRMatrix(const Matrix4d& M) {
  Rcpp::NumericMatrix out(4, 4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      out(i, j) = M(i, j);
  return out;
}

// Regular-step gradient descent on one pyramid level, in the style of ITK's
// RegularStepGradientDescentOptimizerv4: each iteration takes a fixed *physical*
// step down the scale-normalized gradient, and the step is relaxed (halved) only
// when the gradient direction reverses (an overshoot). Stepping is
// unconditional, a momentary uphill move is allowed and self-corrects on the
// next iteration, so the optimizer rides through shallow bumps that the old
// strict reject-and-halve rule would stall on. The best-seen transform is
// tracked and returned, so a level never yields a worse transform than its
// start. The fixed step length (alpha = step / |g|_scaled) bounds each move
// regardless of gradient magnitude, so an uphill step cannot blow up.
void optimizeLevel(LinearProblem& prob, int maxIter,
                   double initStep, double minStep,
                   std::vector<double>& trace,
                   bool verbose = false,
                   const std::string& stageName = "",
                   int levelIdx = 0, int nLevels = 1,
                   double sigma = 0.0, int shrink = 1) {
  if (verbose) {
    Rcpp::Rcout << "[" << stageName << "] level " << (levelIdx+1) << "/" << nLevels
                << " (shrink=" << shrink << ", sigma=" << std::fixed
                << std::setprecision(1) << sigma << ")"
                << ": max " << maxIter << " iterations" << std::endl;
  }

  const int np = prob.xform.nParams();
  VectorXd scales;
  prob.estimateScales(scales);

  const double relax = 0.5;   // step relaxation on gradient reversal (ITK default)

  VectorXd g;
  double cost = prob.evaluate(&g);
  trace.push_back(cost);

  LinearTransform bestX = prob.xform;   // best-seen parameters + their cost
  double bestCost = cost;
  VectorXd gPrev = g;
  double step = initStep;
  int stopCode = 0;   // 0=maxiter, 1=minStep, 2=gradient vanished

  for (int it = 0; it < maxIter; ++it) {
    // scale-weighted gradient magnitude
    double L = 0.0;
    for (int j = 0; j < np; ++j) L += g[j] * g[j] / scales[j];
    L = std::sqrt(L);
    if (L < 1e-15) { stopCode = 2; break; }

    // overshoot detection: reversal of the scaled gradient direction
    if (it > 0) {
      double dot = 0.0;
      for (int j = 0; j < np; ++j) dot += g[j] * gPrev[j] / scales[j];
      if (dot < 0.0) {
        step *= relax;
        if (step < minStep) { stopCode = 1; break; }
      }
    }

    const double alpha = step / L;
    VectorXd delta(np);
    for (int j = 0; j < np; ++j) delta[j] = -alpha * g[j] / scales[j];

    gPrev = g;
    prob.xform.applyStep(delta);          // unconditional regular step
    cost = prob.evaluate(&g);
    trace.push_back(cost);

    if (cost < bestCost) { bestCost = cost; bestX = prob.xform; }

    if (verbose) {
      Rcpp::Rcout << "  it " << std::setw(5) << (it+1)
                  << "  step = " << std::setprecision(6) << step
                  << ": cost = " << std::fixed << std::setprecision(6) << cost
                  << std::endl;
    }

    Rcpp::checkUserInterrupt();
  }

  prob.xform = bestX;   // return the best-seen parameters

  if (verbose) {
    if (stopCode == 1)
      Rcpp::Rcout << "  => converged (step < minStep)  best cost="
                  << std::fixed << std::setprecision(6) << bestCost << std::endl;
    else if (stopCode == 2)
      Rcpp::Rcout << "  => gradient vanished  best cost="
                  << std::fixed << std::setprecision(6) << bestCost << std::endl;
    else
      Rcpp::Rcout << "  => max iterations reached  best cost="
                  << std::fixed << std::setprecision(6) << bestCost << std::endl;
  }
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List register_linear_cpp(const Rcpp::NumericVector& fixed,
                               const Rcpp::IntegerVector& fixedDim,
                               const Rcpp::NumericMatrix& fixedVox2Ras,
                               const Rcpp::NumericVector& moving,
                               const Rcpp::IntegerVector& movingDim,
                               const Rcpp::NumericMatrix& movingVox2Ras,
                               const std::string& type,
                               const std::string& metric,
                               const Rcpp::IntegerVector& shrinkFactors,
                               const Rcpp::NumericVector& smoothingSigmas,
                               const Rcpp::IntegerVector& iterations,
                               const double samplingRate,
                               const double learningRate,
                               const int numberOfBins,
                               const unsigned int seed,
                               const Rcpp::Nullable<Rcpp::NumericMatrix>& initTransform,
                               const Rcpp::Nullable<Rcpp::NumericVector>& fixedMask,
                               const Rcpp::Nullable<Rcpp::NumericVector>& movingMask,
                               const bool verbose = false) {
  const int fnx = fixedDim[0], fny = fixedDim[1], fnz = fixedDim[2];
  const int mnx = movingDim[0], mny = movingDim[1], mnz = movingDim[2];

  const Matrix4d fV2R = asMatrix4(fixedVox2Ras);
  const Matrix4d mV2R = asMatrix4(movingVox2Ras);

  // borrowed full-resolution data
  const double* fixedData = &fixed[0];
  const double* movingData = &moving[0];

  // transform state
  LinearTransform xform;
  xform.mode = (type == "rigid") ? LinearMode::RIGID : LinearMode::AFFINE;
  // center = physical center of the fixed image
  Vector4d cc(((double)fnx - 1.0) / 2.0, ((double)fny - 1.0) / 2.0,
              ((double)fnz - 1.0) / 2.0, 1.0);
  xform.center = (fV2R * cc).head<3>();
  if (initTransform.isNotNull()) {
    Rcpp::NumericMatrix it(initTransform);
    xform.fromRas(asMatrix4(it));
  }

  Metric met = Metric::MEANSQUARES;
  if (metric == "mattes") met = Metric::MATTES;
  else if (metric == "cc") met = Metric::CC;

  // Optional moving mask (full moving grid, reg_real): sampled at mapped points
  // so samples landing outside it are dropped. Built once.
  std::vector<reg_real> movMaskStore;
  RegImage3 movMaskImg;
  const RegImage3* movMaskPtr = nullptr;
  if (movingMask.isNotNull()) {
    Rcpp::NumericVector mm(movingMask);
    movMaskStore.assign(mm.begin(), mm.end());
    movMaskImg = RegImage3(movMaskStore.data(), mnx, mny, mnz, mV2R);
    movMaskPtr = &movMaskImg;
  }
  // Optional fixed mask (full fixed grid, reg_real): shrunk per level below to
  // restrict which fixed voxels are sampled.
  std::vector<reg_real> fixedMaskFull;
  const reg_real* fixedMaskData = nullptr;
  if (fixedMask.isNotNull()) {
    Rcpp::NumericVector fm(fixedMask);
    fixedMaskFull.assign(fm.begin(), fm.end());
    fixedMaskData = fixedMaskFull.data();
  }

  const int nLevels = shrinkFactors.size();
  std::vector<double> fullTrace;

  for (int lev = 0; lev < nLevels; ++lev) {
    const int shrink = shrinkFactors[lev];
    const double sigma = (lev < smoothingSigmas.size()) ? smoothingSigmas[lev] : 0.0;
    const int maxIter = (lev < iterations.size()) ? iterations[lev] : 100;

    // smooth fixed + moving at full resolution (sigma in voxels);
    // gaussianSmooth converts the double R data to reg_real storage
    std::vector<reg_real> sFixed = gaussianSmooth(fixedData, fnx, fny, fnz, sigma);
    std::vector<reg_real> sMoving = gaussianSmooth(movingData, mnx, mny, mnz, sigma);

    FixedLevel flvl = buildFixedLevel(sFixed.data(), fnx, fny, fnz, fV2R, shrink);
    RegImage3 movImg(sMoving.data(), mnx, mny, mnz, mV2R);

    // shrink the fixed mask to this level (trilinear) and threshold to 0/1
    std::vector<reg_real> maskLevel;
    if (fixedMaskData) {
      FixedLevel ml = buildFixedLevel(fixedMaskData, fnx, fny, fnz, fV2R, shrink);
      maskLevel.resize(ml.size());
      for (std::size_t i = 0; i < ml.size(); ++i)
        maskLevel[i] = (ml.data[i] > reg_real(0.5)) ? reg_real(1) : reg_real(0);
    }

    LinearProblem prob;
    prob.fixed = &flvl;
    prob.moving = &movImg;
    prob.movingMask = movMaskPtr;
    prob.fixedMaskLevel = maskLevel.empty() ? nullptr : maskLevel.data();
    prob.xform = xform;
    prob.metric = met;
    prob.numberOfBins = numberOfBins;
    prob.buildSamples(samplingRate, seed + lev);
    prob.prepareRanges();

    const double levelSpacing = flvl.vox2ras.topLeftCorner<3, 3>().determinant();
    const double initStep = (learningRate > 0.0)
        ? learningRate
        : std::cbrt(std::abs(levelSpacing));
    const double minStep = std::max(1e-5, initStep * 1e-3);

    optimizeLevel(prob, maxIter, initStep, minStep, fullTrace,
                  verbose, type, lev, nLevels, sigma, shrink);

    xform = prob.xform;   // carry to next level
  }

  Matrix4d T = xform.toRas();   // fixed RAS -> moving RAS

  return Rcpp::List::create(
    Rcpp::Named("transform") = asRMatrix(T),
    Rcpp::Named("type") = type,
    Rcpp::Named("metric") = metric,
    Rcpp::Named("metric_trace") = Rcpp::wrap(fullTrace)
  );
}
