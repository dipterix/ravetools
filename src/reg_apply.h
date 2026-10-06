#ifndef RAVETOOLS_REG_APPLY_H
#define RAVETOOLS_REG_APPLY_H

// Applying and inverting registration transforms (shared by reg_apply.cpp and
// reg_syn.cpp).
//
// * A composite transform is a chain of stages, each either a 4x4 RAS->RAS
//   affine or a dense RAS displacement field living on its own voxel grid. A
//   point is pushed through the stages in order (stage 1 first), which is the
//   `antsApplyTransforms -t ...` convention: a registration's [warp, affine]
//   list maps a target point r to the source point A(r + u(r)).
// * Fields are sampled trilinearly in their own voxel space and contribute a
//   zero displacement outside their grid (the rule of ITK's displacement-field
//   transform), so a point that leaves the field is still carried by the other
//   stages.
// * invertDisplacementField computes a numerically accurate inverse of a field
//   on its own grid by damped fixed-point iteration (Chen et al. 2008), using
//   the damping schedule and update cap of the ITK inversion filter, and stops
//   on voxel-unit residual tolerances. Every sweep writes each node from one
//   thread and the convergence statistics are reduced serially, so the result
//   is identical for any thread count.
//
// Header-only, in its own namespace, so it cannot collide with the other
// registration translation units.

#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include "reg_interp.h"
#include "TinyParallel.h"

namespace ravereg_apply {

using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// Trilinear sample of a three-component field at continuous 0-based voxel
// coordinates, with the border rule of ITK's displacement-field transform: a
// point is inside the field while it lies within the voxel extent
// [-0.5, n - 0.5) on every axis (in the half-voxel margin around the outer
// nodes the edge value is used, on both the low and the high side), and
// outside that extent (or for non-finite coordinates) the displacement is zero
// and false is returned. Interpolation weights follow the same association
// order as ravereg::trilinearSample, so a sample taken exactly on a node
// returns the node value.
template <typename T>
inline bool sampleField3(const T* ux, const T* uy, const T* uz,
                         const int nx, const int ny, const int nz,
                         const double cx, const double cy, const double cz,
                         double& vx, double& vy, double& vz)
{
  vx = 0.0; vy = 0.0; vz = 0.0;
  if (!(std::isfinite(cx) && std::isfinite(cy) && std::isfinite(cz))) return false;
  if (cx < -0.5 || cy < -0.5 || cz < -0.5 ||
      cx >= static_cast<double>(nx) - 0.5 ||
      cy >= static_cast<double>(ny) - 0.5 ||
      cz >= static_cast<double>(nz) - 0.5) return false;

  // weights from the unclamped lower corner, then each neighbour index is
  // clamped into the grid separately (as ITK's linear interpolator does): in
  // the half-voxel margin on either side both neighbours collapse onto the edge
  // node, so the edge value is returned and nothing is extrapolated
  const int fx0 = static_cast<int>(std::floor(cx));
  const int fy0 = static_cast<int>(std::floor(cy));
  const int fz0 = static_cast<int>(std::floor(cz));
  const double fx = cx - fx0, fy = cy - fy0, fz = cz - fz0;
  auto clampi = [](const int v, const int n) { return (v < 0) ? 0 : ((v > n - 1) ? n - 1 : v); };
  const int x0 = clampi(fx0, nx), x1 = clampi(fx0 + 1, nx);
  const int y0 = clampi(fy0, ny), y1 = clampi(fy0 + 1, ny);
  const int z0 = clampi(fz0, nz), z1 = clampi(fz0 + 1, nz);

  const std::size_t sy = static_cast<std::size_t>(nx);
  const std::size_t sz = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  const std::size_t i000 = static_cast<std::size_t>(x0) + sy * y0 + sz * z0;
  const std::size_t i100 = static_cast<std::size_t>(x1) + sy * y0 + sz * z0;
  const std::size_t i010 = static_cast<std::size_t>(x0) + sy * y1 + sz * z0;
  const std::size_t i110 = static_cast<std::size_t>(x1) + sy * y1 + sz * z0;
  const std::size_t i001 = static_cast<std::size_t>(x0) + sy * y0 + sz * z1;
  const std::size_t i101 = static_cast<std::size_t>(x1) + sy * y0 + sz * z1;
  const std::size_t i011 = static_cast<std::size_t>(x0) + sy * y1 + sz * z1;
  const std::size_t i111 = static_cast<std::size_t>(x1) + sy * y1 + sz * z1;

  auto lerp = [&](const T* d) -> double {
    const double c000 = static_cast<double>(d[i000]), c100 = static_cast<double>(d[i100]);
    const double c010 = static_cast<double>(d[i010]), c110 = static_cast<double>(d[i110]);
    const double c001 = static_cast<double>(d[i001]), c101 = static_cast<double>(d[i101]);
    const double c011 = static_cast<double>(d[i011]), c111 = static_cast<double>(d[i111]);
    const double c00 = c000 + fx * (c100 - c000);
    const double c10 = c010 + fx * (c110 - c010);
    const double c01 = c001 + fx * (c101 - c001);
    const double c11 = c011 + fx * (c111 - c011);
    const double c0 = c00 + fy * (c10 - c00);
    const double c1 = c01 + fy * (c11 - c01);
    return c0 + fz * (c1 - c0);
  };
  vx = lerp(ux); vy = lerp(uy); vz = lerp(uz);
  return true;
}

// One stage of a composite transform (a borrowed view: the caller keeps the
// field memory alive).
struct Stage {
  bool isField = false;
  Matrix4d M = Matrix4d::Identity();          // affine: RAS -> RAS
  const double* ux = nullptr;                 // field components (RAS mm), x fastest
  const double* uy = nullptr;
  const double* uz = nullptr;
  int nx = 0, ny = 0, nz = 0;                 // field grid
  Matrix4d ras2vox = Matrix4d::Identity();    // field grid: RAS -> 0-based voxel
};

// Push a RAS point through the chain (stage 0 first). Returns false when the
// result is not finite (a non-finite input, or a singular stage).
inline bool mapPoint(const Stage* stages, const int nStages, Vector3d& p)
{
  for (int s = 0; s < nStages; ++s) {
    const Stage& st = stages[s];
    const Vector4d h(p[0], p[1], p[2], 1.0);
    if (!st.isField) {
      p = (st.M * h).head<3>();
    } else {
      const Vector4d c = st.ras2vox * h;
      double dx, dy, dz;
      sampleField3(st.ux, st.uy, st.uz, st.nx, st.ny, st.nz, c[0], c[1], c[2], dx, dy, dz);
      p[0] += dx; p[1] += dy; p[2] += dz;
    }
  }
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

// ---------------------------------------------------------------------------
// Fixed-point inversion of a displacement field.
//
// Given u on a grid (RAS displacements), find v on the same grid such that
// (id + u) o (id + v) = id, i.e. v(y) = -u(y + v(y)) at every node y. Starting
// from v = 0, each iteration evaluates the residual e(y) = v(y) + u(y + v(y))
// and moves v against it, v <- v - eps * e, with eps = 0.75 on the first
// iteration and 0.5 afterwards, the per-node step being capped at eps times the
// largest residual (both in voxel units). The iteration stops once the mean and
// the maximum residual over the nodes whose sample point lies inside the grid
// are below the tolerances, or after maxIterations.
// ---------------------------------------------------------------------------
struct InvertOptions {
  int maxIterations = 50;
  double meanTolerance = 1e-5;   // voxel units
  double maxTolerance = 1e-3;    // voxel units
};

struct InvertResult {
  int iterations = 0;
  double meanResidual = 0.0;     // over nodes sampling inside the grid
  double maxResidual = 0.0;
};

// Residual sweep: e(y) = v(y) + u(y + v(y)), plus |e| in voxel units and a
// flag recording whether y + v(y) fell inside the grid.
template <typename T>
struct InvertResidualWorker : public TinyParallel::Worker {
  const T *ux, *uy, *uz;
  int nx, ny, nz;
  Matrix4d vox2ras, ras2vox;
  Matrix3d ras2voxR;
  const double *vx, *vy, *vz;
  double *ex, *ey, *ez, *en;
  unsigned char* inside;

  void operator()(std::size_t begin, std::size_t end) override {
    const std::size_t sz = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    for (std::size_t o = begin; o < end; ++o) {
      const std::size_t k = o / sz;
      const std::size_t r = o - k * sz;
      const std::size_t j = r / static_cast<std::size_t>(nx);
      const std::size_t i = r - j * static_cast<std::size_t>(nx);
      const Vector4d h(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k), 1.0);
      const Vector4d yh = vox2ras * h;
      const Vector4d ph(yh[0] + vx[o], yh[1] + vy[o], yh[2] + vz[o], 1.0);
      const Vector4d c = ras2vox * ph;
      double dx, dy, dz;
      const bool ok = sampleField3(ux, uy, uz, nx, ny, nz, c[0], c[1], c[2], dx, dy, dz);
      const double e0 = vx[o] + dx, e1 = vy[o] + dy, e2 = vz[o] + dz;
      ex[o] = e0; ey[o] = e1; ez[o] = e2;
      en[o] = (ras2voxR * Vector3d(e0, e1, e2)).norm();
      inside[o] = ok ? 1 : 0;
    }
  }
};

// Update sweep: v <- v - eps * e, the voxel-unit step capped at eps * cap.
struct InvertUpdateWorker : public TinyParallel::Worker {
  double *vx, *vy, *vz;
  const double *ex, *ey, *ez, *en;
  double eps, cap;

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      double s = eps;
      if (en[o] > cap && en[o] > 0.0) s *= cap / en[o];
      vx[o] -= s * ex[o];
      vy[o] -= s * ey[o];
      vz[o] -= s * ez[o];
    }
  }
};

// Invert the field (ux, uy, uz) on an nx x ny x nz grid with the given vox2ras
// into (vx, vy, vz), which must each hold nx*ny*nz doubles. Deterministic for
// any thread count.
template <typename T>
inline InvertResult invertDisplacementField(const T* ux, const T* uy, const T* uz,
                                            const int nx, const int ny, const int nz,
                                            const Matrix4d& vox2ras,
                                            double* vx, double* vy, double* vz,
                                            const InvertOptions& opt)
{
  const std::size_t N = static_cast<std::size_t>(nx) *
                        static_cast<std::size_t>(ny) *
                        static_cast<std::size_t>(nz);
  std::fill(vx, vx + N, 0.0);
  std::fill(vy, vy + N, 0.0);
  std::fill(vz, vz + N, 0.0);

  std::vector<double> ex(N), ey(N), ez(N), en(N);
  std::vector<unsigned char> inside(N, 0);

  InvertResidualWorker<T> rw;
  rw.ux = ux; rw.uy = uy; rw.uz = uz;
  rw.nx = nx; rw.ny = ny; rw.nz = nz;
  const Matrix4d ras2vox = vox2ras.inverse();
  rw.vox2ras = vox2ras;
  rw.ras2vox = ras2vox;
  rw.ras2voxR = ras2vox.topLeftCorner<3, 3>();
  rw.vx = vx; rw.vy = vy; rw.vz = vz;
  rw.ex = ex.data(); rw.ey = ey.data(); rw.ez = ez.data(); rw.en = en.data();
  rw.inside = inside.data();

  InvertUpdateWorker uw;
  uw.vx = vx; uw.vy = vy; uw.vz = vz;
  uw.ex = ex.data(); uw.ey = ey.data(); uw.ez = ez.data(); uw.en = en.data();

  InvertResult res;
  const int maxIter = std::max(0, opt.maxIterations);
  for (int it = 0; ; ++it) {
    // main thread only (the sweeps below run the workers): lets a long
    // inversion of a large field be interrupted
    Rcpp::checkUserInterrupt();
    TinyParallel::parallelFor(0, N, rw, 1024);

    // serial, order-independent statistics (keeps the stop decision, hence the
    // field, identical across thread counts)
    double maxAll = 0.0, maxIn = 0.0, sumIn = 0.0;
    std::size_t nIn = 0;
    for (std::size_t o = 0; o < N; ++o) {
      const double e = en[o];
      if (e > maxAll) maxAll = e;
      if (inside[o]) {
        sumIn += e;
        ++nIn;
        if (e > maxIn) maxIn = e;
      }
    }
    res.iterations = it;
    res.meanResidual = (nIn > 0) ? sumIn / static_cast<double>(nIn) : 0.0;
    res.maxResidual = maxIn;
    if (it >= maxIter) break;
    if (maxIn <= opt.maxTolerance && res.meanResidual <= opt.meanTolerance) break;

    const double eps = (it == 0) ? 0.75 : 0.5;
    uw.eps = eps;
    uw.cap = eps * maxAll;
    TinyParallel::parallelFor(0, N, uw, 4096);
  }
  return res;
}

} // namespace ravereg_apply

#endif // RAVETOOLS_REG_APPLY_H
