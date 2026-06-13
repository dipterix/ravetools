#ifndef RAVETOOLS_REG_CORE_H
#define RAVETOOLS_REG_CORE_H

// Core building blocks for native 3D image registration (rigid / affine / SyN).
//
// Everything here works in RAS physical space. Each image carries its own
// vox2ras (0-indexed voxel -> RAS) 4x4 transform, matching the convention used
// by `resample_3d_volume()` / `resample3D()`.
//
// NOTE on handedness: ITK / ANTs operate internally in LPS; this package uses
// RAS throughout. Transforms and displacement fields produced here are RAS->RAS.

#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include "reg_interp.h"   // trilinearSample (Eigen-free, shared with resample3D)
#include "TinyParallel.h"

namespace ravereg {

using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// Bulk voxel storage precision. The geometry (the Eigen 4x4 / 3-vectors above),
// the optimizer parameters, and the metric reductions all stay in double; only
// the large voxel arrays (images, pyramid levels, displacement fields) use this
// reduced precision, roughly halving registration memory. Sampling still reads
// these arrays and accumulates in double, matching ITK's float-image convention.
typedef float reg_real;

// One separable 1D convolution pass along `axis`, parallelized over the
// orthogonal "lines". Reads `in`, writes `out` (distinct buffers). Borders
// clamp to the edge value. Values are reg_real; the accumulator is double.
struct SepConvWorker : public TinyParallel::Worker {
  int axis, nx, ny, nz, radius;
  const double* kernel;
  const reg_real* in;
  reg_real* out;

  SepConvWorker(int axis_, int nx_, int ny_, int nz_, int radius_,
                const double* kernel_, const reg_real* in_, reg_real* out_)
    : axis(axis_), nx(nx_), ny(ny_), nz(nz_), radius(radius_),
      kernel(kernel_), in(in_), out(out_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t sy = static_cast<std::size_t>(nx);
    const std::size_t sz = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    for (std::size_t L = begin; L < end; ++L) {
      int a, b, n, stride;
      std::size_t base;
      if (axis == 0) {            // x lines, indexed by (y, z)
        a = (int)(L % ny); b = (int)(L / ny);
        base = sy * a + sz * b; n = nx; stride = 1;
      } else if (axis == 1) {     // y lines, indexed by (x, z)
        a = (int)(L % nx); b = (int)(L / nx);
        base = a + sz * b; n = ny; stride = (int)sy;
      } else {                    // z lines, indexed by (x, y)
        a = (int)(L % nx); b = (int)(L / nx);
        base = a + sy * b; n = nz; stride = (int)sz;
      }
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (int t = -radius; t <= radius; ++t) {
          int ii = i + t; if (ii < 0) ii = 0; else if (ii >= n) ii = n - 1;
          acc += kernel[t + radius] * in[base + (std::size_t)stride * ii];
        }
        out[base + (std::size_t)stride * i] = acc;
      }
    }
  }
};

// ---------------------------------------------------------------------------
// Image view: a borrowed 3D volume plus its vox2ras / ras2vox transforms
// ---------------------------------------------------------------------------
struct RegImage3 {
  const reg_real* data = nullptr;
  int nx = 0, ny = 0, nz = 0;
  Matrix4d vox2ras = Matrix4d::Identity();
  Matrix4d ras2vox = Matrix4d::Identity();
  Matrix3d ras2voxR = Matrix3d::Identity();   // upper-left 3x3 of ras2vox

  RegImage3() = default;

  RegImage3(const reg_real* d, int x, int y, int z, const Matrix4d& v2r)
    : data(d), nx(x), ny(y), nz(z) {
    setVox2Ras(v2r);
  }

  void setVox2Ras(const Matrix4d& v2r) {
    vox2ras = v2r;
    ras2vox = v2r.inverse();
    ras2voxR = ras2vox.topLeftCorner<3, 3>();
  }

  std::size_t size() const {
    return static_cast<std::size_t>(nx) *
           static_cast<std::size_t>(ny) *
           static_cast<std::size_t>(nz);
  }

  // RAS coordinate of a 0-based voxel index.
  inline Vector3d voxelToRas(double i, double j, double k) const {
    Vector4d v(i, j, k, 1.0);
    return (vox2ras * v).head<3>();
  }

  // Sample the volume at a RAS point. Optionally returns the spatial gradient
  // d value / d RAS (3 components). Returns false if the point is out of bounds.
  inline bool sampleRas(const Vector3d& p, double& value, Vector3d* gradRas = nullptr) const {
    Vector4d ph(p[0], p[1], p[2], 1.0);
    Vector4d v = ras2vox * ph;
    if (gradRas) {
      double gvox[3];
      bool ok = trilinearSample(data, nx, ny, nz, v[0], v[1], v[2], value, gvox);
      if (ok) {
        // d value / d RAS = ras2voxR^T * (d value / d voxel)
        Vector3d gv(gvox[0], gvox[1], gvox[2]);
        *gradRas = ras2voxR.transpose() * gv;
      } else {
        gradRas->setZero();
      }
      return ok;
    }
    return trilinearSample(data, nx, ny, nz, v[0], v[1], v[2], value, nullptr);
  }

  // Mean voxel spacing (geometric), useful as a length scale per level.
  double meanSpacing() const {
    double det = std::abs(vox2ras.topLeftCorner<3, 3>().determinant());
    return std::cbrt(det);
  }
};

// ---------------------------------------------------------------------------
// Separable Gaussian smoothing (sigma in voxels). Borders are clamped.
//
// gaussianSmoothInto writes the result into caller-provided `dst`, using
// `scratch` (both size nx*ny*nz) as the ping-pong buffer -- no allocation, so
// hot per-iteration callers (the SyN field smoothing) can reuse buffers. The
// three separable passes match the allocating wrapper exactly, so the output is
// bit-identical. `src == dst` is allowed (src is read only into scratch first).
// ---------------------------------------------------------------------------
template <typename T>
inline void gaussianSmoothInto(const T* src, reg_real* dst, reg_real* scratch,
                               int nx, int ny, int nz, double sigma)
{
  const std::size_t n = static_cast<std::size_t>(nx) *
                        static_cast<std::size_t>(ny) *
                        static_cast<std::size_t>(nz);
  if (sigma <= 0.0) {
    for (std::size_t i = 0; i < n; ++i) dst[i] = static_cast<reg_real>(src[i]);
    return;
  }

  const int radius = std::max(1, static_cast<int>(std::ceil(3.0 * sigma)));
  std::vector<double> k(2 * radius + 1);
  const double s2 = 2.0 * sigma * sigma;
  double ksum = 0.0;
  for (int i = -radius; i <= radius; ++i) {
    const double w = std::exp(-static_cast<double>(i * i) / s2);
    k[i + radius] = w;
    ksum += w;
  }
  for (double& w : k) w /= ksum;

  // seed scratch with src, then ping-pong scratch <-> dst; result ends in dst
  // (passes: scratch->dst, dst->scratch, scratch->dst).
  for (std::size_t i = 0; i < n; ++i) scratch[i] = static_cast<reg_real>(src[i]);

  SepConvWorker wx{0, nx, ny, nz, radius, k.data(), scratch, dst};
  TinyParallel::parallelFor(0, static_cast<std::size_t>(ny) * nz, wx, 32);
  SepConvWorker wy{1, nx, ny, nz, radius, k.data(), dst, scratch};
  TinyParallel::parallelFor(0, static_cast<std::size_t>(nx) * nz, wy, 32);
  SepConvWorker wz{2, nx, ny, nz, radius, k.data(), scratch, dst};
  TinyParallel::parallelFor(0, static_cast<std::size_t>(nx) * ny, wz, 32);
}

// Allocating convenience wrapper (used by the non-hot per-level image
// smoothing). Identical output to gaussianSmoothInto.
template <typename T>
inline std::vector<reg_real> gaussianSmooth(const T* data,
                                            int nx, int ny, int nz,
                                            double sigma)
{
  const std::size_t n = static_cast<std::size_t>(nx) *
                        static_cast<std::size_t>(ny) *
                        static_cast<std::size_t>(nz);
  std::vector<reg_real> dst(n), scratch(n);
  gaussianSmoothInto(data, dst.data(), scratch.data(), nx, ny, nz, sigma);
  return dst;
}

// ---------------------------------------------------------------------------
// Multiresolution pyramid level for the FIXED image.
//
// Holds a shrunk + smoothed copy of the fixed volume, with its own vox2ras so
// that each coarse voxel maps to the correct RAS location. Coarse index j maps
// to fine index (shrink * j + (shrink - 1) / 2) -- cell-centered, matching the
// behavior of ITK's ShrinkImageFilter.
// ---------------------------------------------------------------------------
struct FixedLevel {
  std::vector<reg_real> data;
  int nx = 0, ny = 0, nz = 0;
  Matrix4d vox2ras = Matrix4d::Identity();

  std::size_t size() const {
    return static_cast<std::size_t>(nx) *
           static_cast<std::size_t>(ny) *
           static_cast<std::size_t>(nz);
  }
};

// Build a fixed-image pyramid level from a (already smoothed) full-resolution
// fixed volume and the full-resolution vox2ras.
inline FixedLevel buildFixedLevel(const reg_real* smoothFixed,
                                  int fnx, int fny, int fnz,
                                  const Matrix4d& fixedVox2Ras,
                                  int shrink)
{
  FixedLevel lvl;
  if (shrink < 1) shrink = 1;
  lvl.nx = std::max(1, static_cast<int>(std::round(static_cast<double>(fnx) / shrink)));
  lvl.ny = std::max(1, static_cast<int>(std::round(static_cast<double>(fny) / shrink)));
  lvl.nz = std::max(1, static_cast<int>(std::round(static_cast<double>(fnz) / shrink)));
  lvl.data.assign(lvl.size(), 0.0);

  const double off = (shrink - 1) / 2.0;

  // coarse vox2ras = fixedVox2Ras * S, where S maps coarse idx -> fine idx
  Matrix4d S = Matrix4d::Identity();
  S(0, 0) = shrink; S(1, 1) = shrink; S(2, 2) = shrink;
  S(0, 3) = off;    S(1, 3) = off;    S(2, 3) = off;
  lvl.vox2ras = fixedVox2Ras * S;

  const std::size_t sy = static_cast<std::size_t>(lvl.nx);
  const std::size_t sz = static_cast<std::size_t>(lvl.nx) * static_cast<std::size_t>(lvl.ny);
  for (int k = 0; k < lvl.nz; ++k) {
    const double fk = shrink * k + off;
    for (int j = 0; j < lvl.ny; ++j) {
      const double fj = shrink * j + off;
      for (int i = 0; i < lvl.nx; ++i) {
        const double fi = shrink * i + off;
        double val = 0.0;
        trilinearSample(smoothFixed, fnx, fny, fnz, fi, fj, fk, val, nullptr);
        lvl.data[i + sy * j + sz * k] = val;
      }
    }
  }
  return lvl;
}

} // namespace ravereg

#endif // RAVETOOLS_REG_CORE_H
