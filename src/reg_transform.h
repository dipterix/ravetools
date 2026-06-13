#ifndef RAVETOOLS_REG_TRANSFORM_H
#define RAVETOOLS_REG_TRANSFORM_H

// Linear transform parameterization for registration.
//
// The transform maps a FIXED-image RAS point p to a MOVING-image RAS point:
//
//     q = A * (p - center) + center + t
//
// `center` is fixed (the physical center of the fixed image), so rotations act
// about the image center. The optimizer works with incremental parameters:
//   - AFFINE: 12 params = matrix entries (row-major, 9) then translation (3),
//     applied additively. Matches ITK's AffineTransform parameter ordering.
//   - RIGID : 6 params = rotation vector in so(3) (3) then translation (3).
//     The rotation increment is applied on the left via the matrix exponential
//     (Rodrigues), A <- exp([w]_x) * A; translation is additive.

#include <RcppEigen.h>
#include <cmath>

namespace ravereg {

using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;
using Eigen::VectorXd;

enum class LinearMode { RIGID, AFFINE };

// Rodrigues: rotation matrix for rotation vector w (axis * angle).
inline Matrix3d expSO3(const Vector3d& w) {
  const double theta = w.norm();
  Matrix3d K = Matrix3d::Zero();
  if (theta < 1e-12) {
    // first-order: I + [w]_x
    K(0, 1) = -w[2]; K(0, 2) =  w[1];
    K(1, 0) =  w[2]; K(1, 2) = -w[0];
    K(2, 0) = -w[1]; K(2, 1) =  w[0];
    return Matrix3d::Identity() + K;
  }
  const Vector3d a = w / theta;
  K(0, 1) = -a[2]; K(0, 2) =  a[1];
  K(1, 0) =  a[2]; K(1, 2) = -a[0];
  K(2, 0) = -a[1]; K(2, 1) =  a[0];
  return Matrix3d::Identity() + std::sin(theta) * K +
         (1.0 - std::cos(theta)) * (K * K);
}

struct LinearTransform {
  LinearMode mode = LinearMode::AFFINE;
  Matrix3d A = Matrix3d::Identity();
  Vector3d t = Vector3d::Zero();
  Vector3d center = Vector3d::Zero();

  int nParams() const { return mode == LinearMode::RIGID ? 6 : 12; }

  inline Vector3d map(const Vector3d& p) const {
    return A * (p - center) + center + t;
  }

  // 4x4 RAS->RAS transform mapping fixed RAS -> moving RAS.
  Matrix4d toRas() const {
    Matrix4d M = Matrix4d::Identity();
    M.topLeftCorner<3, 3>() = A;
    M.topRightCorner<3, 1>() = center - A * center + t;
    return M;
  }

  // Initialize A, t from a 4x4 RAS->RAS transform (center must already be set).
  void fromRas(const Matrix4d& M) {
    A = M.topLeftCorner<3, 3>();
    Vector3d offset = M.topRightCorner<3, 1>();
    // offset = center - A*center + t  =>  t = offset - center + A*center
    t = offset - center + A * center;
  }

  // Apply an incremental parameter step (length nParams()).
  void applyStep(const VectorXd& d) {
    if (mode == LinearMode::RIGID) {
      Vector3d w(d[0], d[1], d[2]);
      A = expSO3(w) * A;
      t += Vector3d(d[3], d[4], d[5]);
    } else {
      // additive on matrix (row-major) then translation
      A(0, 0) += d[0]; A(0, 1) += d[1]; A(0, 2) += d[2];
      A(1, 0) += d[3]; A(1, 1) += d[4]; A(1, 2) += d[5];
      A(2, 0) += d[6]; A(2, 1) += d[7]; A(2, 2) += d[8];
      t += Vector3d(d[9], d[10], d[11]);
    }
  }
};

} // namespace ravereg

#endif // RAVETOOLS_REG_TRANSFORM_H
