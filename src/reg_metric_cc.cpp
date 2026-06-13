// Global normalized cross-correlation metric for linear registration.
//
// This is ITK's CorrelationImageToImageMetricv4 (a.k.a. ANTs "GC"): the
// normalized correlation taken over all sampled voxels. It is well suited to
// same-modality linear alignment and has a clean analytic gradient. (The local
// neighborhood CC used by ANTs SyN is implemented separately for the deformable
// phase.)
//
//   NCC = Sfm / sqrt(Sff * Smm),   a_i = F_i - meanF,  b_i = M_i - meanM
//   Sfm = sum a_i b_i,  Sff = sum a_i^2,  Smm = sum b_i^2
//   cost = -NCC
//
//   d(NCC)/dmu = (1 / sqrt(Sff*Smm)) * sum_i [ a_i - (Sfm/Smm) b_i ] * dM_i/dmu

#include "reg_metric.h"

namespace ravereg {

double LinearProblem::evalCC(VectorXd* grad) {
  const int np = xform.nParams();
  if (grad) grad->setZero(np);

  const std::size_t ns = sampleIdx.size();

  double sumF = 0.0, sumM = 0.0;
  long n = 0;
  for (std::size_t s = 0; s < ns; ++s) {
    if (!sValid[s]) continue;
    sumF += sF[s];
    sumM += sMw[s];
    ++n;
  }
  if (n < 2) return 0.0;

  const double meanF = sumF / n;
  const double meanM = sumM / n;

  double Sff = 0.0, Smm = 0.0, Sfm = 0.0;
  for (std::size_t s = 0; s < ns; ++s) {
    if (!sValid[s]) continue;
    const double a = sF[s] - meanF;
    const double b = sMw[s] - meanM;
    Sff += a * a;
    Smm += b * b;
    Sfm += a * b;
  }
  const double denom = std::sqrt(Sff * Smm);
  if (denom < 1e-12) return 0.0;
  const double ncc = Sfm / denom;

  if (grad) {
    const double invDenom = 1.0 / denom;
    const double sfmOverSmm = Sfm / Smm;
    for (std::size_t s = 0; s < ns; ++s) {
      if (!sValid[s]) continue;
      const Vector3d yv(sPx[s] - xform.center[0],
                        sPy[s] - xform.center[1],
                        sPz[s] - xform.center[2]);
      const Vector3d gM(sGx[s], sGy[s], sGz[s]);
      const double a = sF[s] - meanF;
      const double b = sMw[s] - meanM;
      const double coef = invDenom * (a - sfmOverSmm * b);
      // cost = -NCC, so residual factor is -coef
      accumulateLinearGrad(*grad, -coef, yv, gM);
    }
  }

  return -ncc;
}

} // namespace ravereg
