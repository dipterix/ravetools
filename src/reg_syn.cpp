// Deformable registration: greedy symmetric diffeomorphic ("SyN"-style).
//
// After a linear (affine) stage, a dense RAS displacement field on the FIXED
// grid captures the residual nonlinear deformation. Each iteration:
//   1. warp the moving image to the fixed grid through affine o (id + disp);
//   2. form a similarity force (local normalized CC, or mean squares);
//   3. Gaussian-smooth the update (flow_sigma), the fluid-like regularizer;
//   4. take a capped step (grad_step) on the forward field;
//   5. optionally Gaussian-smooth the total field (total_sigma), elastic;
//   6. optionally end the level once the cost of the first channel stops
//      decreasing (the convergence rule of ANTs, ConvergenceMonitor below).
// After the last level the inverse field is computed from the final forward
// field on the full fixed grid by fixed-point inversion (reg_apply.h), so
// forward o inverse and inverse o forward are both identity to within the
// inversion tolerance (the negated forward field used before was only a
// first-order inverse).
//
// Everything is in RAS. This is the functional-equivalence counterpart to ANTs
// SyN; it is diffeomorphic for smooth/moderate deformations rather than a full
// time-integrated velocity field.
//
// Multi-channel (multivariate): several already co-registered modalities can
// drive the same field. Each channel k carries its own fixed/moving pair and
// metric; the per-channel similarity forces are combined as a weighted sum
//   force = sum_k weight[k] * force_k,        (weights normalized to sum 1)
// matching ANTs' multivariate SyN. The reported metric is the same weighted sum
// of the per-channel costs. Channels share the grid, the affine, and the field;
// only the image value/gradient differ per channel.
//
// Performance: every per-iteration sweep over the level grid is parallelized
// with TinyParallel (the warp, the separable box-sums, the CC / mean-squares
// force, the step-cap reduction, the field update), mirroring the parallel
// sampling in the linear stage (reg_metric.h MovingSampler) and the separable
// Gaussian (reg_core.h SepConvWorker). The deformable stage processes every
// voxel each iteration (no metric subsampling, a dense force field and the CC
// box-sums need the full warped image), so parallelism is the main lever. The
// per-voxel force `factor` has no cross-voxel dependence, so the displacement
// field is bit-identical regardless of thread count. The scalar cost is summed
// over fixed blocks that are combined in block order, so the cost trace, and
// with it the convergence stop, is identical for any thread count as well (for
// a given starting transform: the Mattes linear stages still sum in a
// thread-dependent order).

#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <random>
#include <cstdint>
#include "reg_core.h"
#include "reg_apply.h"   // ravereg_apply::invertDisplacementField

using namespace ravereg;

namespace {

struct Field {                 // displacement field on a fixed grid (RAS)
  int nx = 0, ny = 0, nz = 0;
  Matrix4d vox2ras = Matrix4d::Identity();
  std::vector<reg_real> x, y, z;
  std::size_t size() const {
    return (std::size_t)nx * (std::size_t)ny * (std::size_t)nz;
  }
  void alloc(int a, int b, int c) {
    nx = a; ny = b; nz = c;
    x.assign(size(), 0.0); y.assign(size(), 0.0); z.assign(size(), 0.0);
  }
  void zero() {
    std::fill(x.begin(), x.end(), (reg_real)0);
    std::fill(y.begin(), y.end(), (reg_real)0);
    std::fill(z.begin(), z.end(), (reg_real)0);
  }
};

Matrix4d asMatrix4(const Rcpp::NumericMatrix& m) {
  Matrix4d M = Matrix4d::Identity();
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) M(i, j) = m(i, j);
  return M;
}

// Decode a column-major (x fastest) flat index into (i, j, k).
inline void decodeIndex(std::size_t o, int nx, std::size_t sz,
                        int& i, int& j, int& k) {
  const std::size_t kk = o / sz;
  const std::size_t r  = o - kk * sz;
  const std::size_t jj = r / (std::size_t)nx;
  const std::size_t ii = r - jj * (std::size_t)nx;
  i = (int)ii; j = (int)jj; k = (int)kk;
}

// ---------------------------------------------------------------------------
// Separable box-sum over a (2r+1)^3 window with clamped borders, parallelized
// per orthogonal line (one of the three sweeps below runs per axis, with a
// barrier between sweeps). Accepts reg_real or double input; always accumulates
// and returns double so the CC variance terms (which subtract nearly equal
// large sums) stay accurate. Mirrors reg_core.h's SepConvWorker.
// ---------------------------------------------------------------------------
template <typename T>
struct BoxSum1D : public TinyParallel::Worker {
  int axis, nx, ny, nz, r;
  const T* in;
  double* out;
  BoxSum1D(int axis_, int nx_, int ny_, int nz_, int r_,
           const T* in_, double* out_)
    : axis(axis_), nx(nx_), ny(ny_), nz(nz_), r(r_), in(in_), out(out_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    const std::size_t sy = (std::size_t)nx;
    const std::size_t sz = (std::size_t)nx * (std::size_t)ny;
    for (std::size_t L = begin; L < end; ++L) {
      int n, stride; std::size_t base;
      if (axis == 0) {            // x lines, indexed by (y, z)
        int a = (int)(L % ny), b = (int)(L / ny);
        base = sy * a + sz * b; n = nx; stride = 1;
      } else if (axis == 1) {     // y lines, indexed by (x, z)
        int a = (int)(L % nx), b = (int)(L / nx);
        base = a + sz * b; n = ny; stride = (int)sy;
      } else {                    // z lines, indexed by (x, y)
        int a = (int)(L % nx), b = (int)(L / nx);
        base = a + sy * b; n = nz; stride = (int)sz;
      }
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (int t = -r; t <= r; ++t) {
          int ii = i + t; if (ii < 0) ii = 0; else if (ii >= n) ii = n - 1;
          acc += (double)in[base + (std::size_t)stride * ii];
        }
        out[base + (std::size_t)stride * i] = acc;
      }
    }
  }
};

// Box-sum writing into caller-provided `out` (size N), using `scratch` (size N)
// as the intermediate buffer, no allocation, so the per-iteration CC sums
// reuse hoisted buffers. Result lands in `out`.
template <typename T>
void boxSumInto(const T* in, int nx, int ny, int nz, int r,
                std::vector<double>& out, std::vector<double>& scratch) {
  BoxSum1D<T> wx(0, nx, ny, nz, r, in, out.data());
  TinyParallel::parallelFor(0, (std::size_t)ny * nz, wx, 32);
  BoxSum1D<double> wy(1, nx, ny, nz, r, out.data(), scratch.data());
  TinyParallel::parallelFor(0, (std::size_t)nx * nz, wy, 32);
  BoxSum1D<double> wz(2, nx, ny, nz, r, scratch.data(), out.data());
  TinyParallel::parallelFor(0, (std::size_t)nx * ny, wz, 32);
}

// Allocating convenience wrapper (used for the per-level fixed sums, which are
// computed once per level rather than per iteration).
template <typename T>
std::vector<double> boxSum(const std::vector<T>& in, int nx, int ny, int nz, int r) {
  std::vector<double> out(in.size()), scratch(in.size());
  boxSumInto(in.data(), nx, ny, nz, r, out, scratch);
  return out;
}

// ---------------------------------------------------------------------------
// Resample a displacement field onto a new fixed grid (trilinear, per RAS
// component; displacements are in RAS so they transfer without rescaling).
// Parallelized over the destination grid (each output voxel independent).
// ---------------------------------------------------------------------------
struct ResampleFieldWorker : public TinyParallel::Worker {
  Matrix4d dstV2R, s2v;
  const reg_real *sx, *sy, *sz;
  int snx, sny, snz;
  int dnx; std::size_t szd;
  reg_real *dx, *dy, *dz;

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      int i, j, k; decodeIndex(o, dnx, szd, i, j, k);
      Vector4d ph(i, j, k, 1.0);
      Vector3d ras = (dstV2R * ph).head<3>();
      Vector4d rh(ras[0], ras[1], ras[2], 1.0);
      Vector4d sv = s2v * rh;
      double vx = 0, vy = 0, vz = 0;
      trilinearSample(sx, snx, sny, snz, sv[0], sv[1], sv[2], vx);
      trilinearSample(sy, snx, sny, snz, sv[0], sv[1], sv[2], vy);
      trilinearSample(sz, snx, sny, snz, sv[0], sv[1], sv[2], vz);
      dx[o] = vx; dy[o] = vy; dz[o] = vz;
    }
  }
};

void resampleField(const Field& src, Field& dst) {
  ResampleFieldWorker w;
  w.dstV2R = dst.vox2ras;
  w.s2v = src.vox2ras.inverse();
  w.sx = src.x.data(); w.sy = src.y.data(); w.sz = src.z.data();
  w.snx = src.nx; w.sny = src.ny; w.snz = src.nz;
  w.dnx = dst.nx; w.szd = (std::size_t)dst.nx * (std::size_t)dst.ny;
  w.dx = dst.x.data(); w.dy = dst.y.data(); w.dz = dst.z.data();
  TinyParallel::parallelFor(0, dst.size(), w, 256);
}

// Smooth each field component in place through a caller-provided scratch buffer
// (size = field size), so the hot per-iteration smoothing allocates nothing.
// Bit-identical to allocating the result via gaussianSmooth.
void smoothFieldInto(Field& f, double sigma, reg_real* scratch) {
  if (sigma <= 0.0) return;
  gaussianSmoothInto(f.x.data(), f.x.data(), scratch, f.nx, f.ny, f.nz, sigma);
  gaussianSmoothInto(f.y.data(), f.y.data(), scratch, f.nx, f.ny, f.nz, sigma);
  gaussianSmoothInto(f.z.data(), f.z.data(), scratch, f.nx, f.ny, f.nz, sigma);
}

// ---------------------------------------------------------------------------
// Per-level base RAS grid: x = levelV2R * (i, j, k). Constant across all
// iterations of a level, so it is precomputed once and reused by the warp.
// ---------------------------------------------------------------------------
struct BaseGridWorker : public TinyParallel::Worker {
  Matrix4d v2r;
  int nx; std::size_t sz;
  double *bx, *by, *bz;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      int i, j, k; decodeIndex(o, nx, sz, i, j, k);
      Vector4d ph(i, j, k, 1.0);
      Vector3d x = (v2r * ph).head<3>();
      bx[o] = x[0]; by[o] = x[1]; bz[o] = x[2];
    }
  }
};

// ---------------------------------------------------------------------------
// Warp one moving channel onto the level grid through aff o (id + disp), and
// record the warped value plus the fixed-space spatial gradient gF = affR^T gM.
// Each output index is written by exactly one thread (no reduction). Mirrors
// reg_metric.h's MovingSampler.
// ---------------------------------------------------------------------------
struct WarpWorker : public TinyParallel::Worker {
  const RegImage3* mov;
  const RegImage3* movMask = nullptr;   // optional moving mask
  const reg_real* warpMask = nullptr;   // optional fixed mask, dilated by ccRadius
  Matrix4d aff;
  Matrix3d affRT;
  const double *bx, *by, *bz;          // base RAS grid (constant per level)
  const reg_real *dispx, *dispy, *dispz;  // current forward displacement
  double *Mw, *gFx, *gFy, *gFz;
  // optional: 1 where the sample lies inside the moving image and mask (the
  // samples of the mutual information, and in the symmetric path the samples
  // that drive either half), 0 elsewhere
  uint8_t* valid = nullptr;

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      // Outside the dilated fixed mask there is nothing to drive and no force
      // voxel needs this value in its window: skip the (dominant) sampling work.
      if (warpMask && warpMask[o] == reg_real(0)) {
        Mw[o] = 0.0; gFx[o] = 0.0; gFy[o] = 0.0; gFz[o] = 0.0;
        if (valid) valid[o] = 0;
        continue;
      }
      Vector4d xdh(bx[o] + dispx[o], by[o] + dispy[o], bz[o] + dispz[o], 1.0);
      Vector3d q = (aff * xdh).head<3>();
      double mval; Vector3d gM;
      const bool inside = mov->sampleRas(q, mval, &gM);
      // A point mapping outside the moving mask is treated as background (no
      // value, no gradient) so it contributes neither force nor structure.
      if (movMask) {
        double mv; movMask->sampleRas(q, mv);
        if (mv <= 0.5) {
          Mw[o] = 0.0; gFx[o] = 0.0; gFy[o] = 0.0; gFz[o] = 0.0;
          if (valid) valid[o] = 0;
          continue;
        }
      }
      Mw[o] = mval;
      Vector3d gF = affRT * gM;
      gFx[o] = gF[0]; gFy[o] = gF[1]; gFz[o] = gF[2];
      if (valid) valid[o] = inside ? 1 : 0;
    }
  }
};

// Fused local-CC windowed sums. Local CC needs box-sums of M, M^2 and F*M over
// the (2r+1)^3 window each iteration (F, F^2 are constant and precomputed). The
// x-pass forms all three running sums straight from Mw and F in one traversal
// (no MM/FM materialization); the y/z passes carry the three partial sums
// together. This reads Mw/F once instead of three times and cuts the launches
// from ~11 (three boxSumInto + a product pass) to 3. Output is bit-identical to
// the separate box-sums: each element is the same double sum in the same order.

// x-pass: windowed sums of M, M^2, F*M along x (lines indexed by (y, z)).
struct CCSumsXWorker : public TinyParallel::Worker {
  const double* Mw;
  const reg_real* F;
  int nx, ny, nz, r;
  double *oM, *oMM, *oFM;
  void operator()(std::size_t begin, std::size_t end) override {
    const std::size_t sy = (std::size_t)nx, sz = (std::size_t)nx * ny;
    for (std::size_t L = begin; L < end; ++L) {
      const int a = (int)(L % ny), b = (int)(L / ny);
      const std::size_t base = sy * a + sz * b;
      for (int x = 0; x < nx; ++x) {
        double sM = 0.0, sMM = 0.0, sFM = 0.0;
        for (int t = -r; t <= r; ++t) {
          int xx = x + t; if (xx < 0) xx = 0; else if (xx >= nx) xx = nx - 1;
          const std::size_t o = base + xx;
          const double m = Mw[o];
          // round the products to double *before* accumulating, matching the old
          // store-to-MM/FM-then-box-sum path; this also stops the compiler from
          // contracting `sum += a*b` into a single-rounding FMA, so the result is
          // bit-identical to the pre-fusion code.
          const double mm = m * m;
          const double fm = (double)F[o] * m;
          sM += m; sMM += mm; sFM += fm;
        }
        oM[base + x] = sM; oMM[base + x] = sMM; oFM[base + x] = sFM;
      }
    }
  }
};

// y/z-pass: windowed sum of three double arrays together along axis (1=y, 2=z).
struct Box3Worker : public TinyParallel::Worker {
  int axis, nx, ny, nz, r;
  const double *iA, *iB, *iC;
  double *oA, *oB, *oC;
  void operator()(std::size_t begin, std::size_t end) override {
    const std::size_t sy = (std::size_t)nx, sz = (std::size_t)nx * ny;
    for (std::size_t L = begin; L < end; ++L) {
      int n, stride; std::size_t base;
      if (axis == 1) {            // y lines, indexed by (x, z)
        const int a = (int)(L % nx), b = (int)(L / nx);
        base = a + sz * b; n = ny; stride = (int)sy;
      } else {                    // z lines, indexed by (x, y)
        const int a = (int)(L % nx), b = (int)(L / nx);
        base = a + sy * b; n = nz; stride = (int)sz;
      }
      for (int i = 0; i < n; ++i) {
        double sA = 0.0, sB = 0.0, sC = 0.0;
        for (int t = -r; t <= r; ++t) {
          int ii = i + t; if (ii < 0) ii = 0; else if (ii >= n) ii = n - 1;
          const std::size_t o = base + (std::size_t)stride * ii;
          sA += iA[o]; sB += iB[o]; sC += iC[o];
        }
        const std::size_t oo = base + (std::size_t)stride * i;
        oA[oo] = sA; oB[oo] = sB; oC[oo] = sC;
      }
    }
  }
};

// Drive the three fused passes; results land in sumM/sumMM/sumFM (scr* are the
// hoisted intermediates). Mirrors boxSumInto's ping-pong (x: in->a, y: a->b,
// z: b->a) but for the three CC quantities at once.
void fusedLocalCCSums(const double* Mw, const reg_real* F,
                      int nx, int ny, int nz, int r,
                      std::vector<double>& sumM, std::vector<double>& sumMM,
                      std::vector<double>& sumFM, std::vector<double>& scrM,
                      std::vector<double>& scrMM, std::vector<double>& scrFM) {
  CCSumsXWorker wx;
  wx.Mw = Mw; wx.F = F; wx.nx = nx; wx.ny = ny; wx.nz = nz; wx.r = r;
  wx.oM = sumM.data(); wx.oMM = sumMM.data(); wx.oFM = sumFM.data();
  TinyParallel::parallelFor(0, (std::size_t)ny * nz, wx, 32);

  Box3Worker wy;
  wy.axis = 1; wy.nx = nx; wy.ny = ny; wy.nz = nz; wy.r = r;
  wy.iA = sumM.data(); wy.iB = sumMM.data(); wy.iC = sumFM.data();
  wy.oA = scrM.data(); wy.oB = scrMM.data(); wy.oC = scrFM.data();
  TinyParallel::parallelFor(0, (std::size_t)nx * nz, wy, 32);

  Box3Worker wz;
  wz.axis = 2; wz.nx = nx; wz.ny = ny; wz.nz = nz; wz.r = r;
  wz.iA = scrM.data(); wz.iB = scrMM.data(); wz.iC = scrFM.data();
  wz.oA = sumM.data(); wz.oB = sumMM.data(); wz.oC = sumFM.data();
  TinyParallel::parallelFor(0, (std::size_t)nx * ny, wz, 32);
}

// Voxels per cost-reduction block. Fixed, so the block partial sums of the cost
// (hence the cost trace and the convergence stop) do not depend on how
// TinyParallel splits the block range among threads (the convention of
// imaging_segment_gmm.cpp).
const std::size_t kCostBlock = 2048;

// Local normalized-CC force, accumulated (weighted) into the shared force
// field; the local-CC sum and the number of valid windows are reduced per block
// of kCostBlock voxels for the cost trace. Each index i is written by one
// thread (disjoint blocks), so the force is bit-identical to the serial result.
struct CCForceWorker : public TinyParallel::Worker {
  const reg_real* F = nullptr;
  const double *Mw = nullptr, *sumF = nullptr, *sumFF = nullptr, *sumM = nullptr,
               *sumMM = nullptr, *sumFM = nullptr;
  const double *gFx = nullptr, *gFy = nullptr, *gFz = nullptr;
  reg_real *fx = nullptr, *fy = nullptr, *fz = nullptr;
  double winN = 0.0, tau = 0.0, wk = 0.0;
  const reg_real* mask = nullptr;   // optional fixed mask: skip force where 0
  std::size_t N = 0;
  double* blockAcc = nullptr;       // per block: sum of the local CC^2
  long* blockCnt = nullptr;         // per block: number of valid windows

  void operator()(std::size_t bBegin, std::size_t bEnd) override {
    for (std::size_t b = bBegin; b < bEnd; ++b) {
      const std::size_t begin = b * kCostBlock;
      const std::size_t end = std::min(N, begin + kCostBlock);
      double accCC = 0.0; long cnt = 0;
      for (std::size_t i = begin; i < end; ++i) {
        if (mask && mask[i] == reg_real(0)) continue;   // outside fixed mask
        const double mF = sumF[i] / winN, mM = sumM[i] / winN;
        const double Sff = sumFF[i] - winN * mF * mF;
        const double Smm = sumMM[i] - winN * mM * mM;
        const double Sfm = sumFM[i] - winN * mF * mM;
        const double denom = Sff * Smm;
        // require both windows to carry real structure (skip flat/background)
        if (Sff > tau && Smm > tau && denom > 1e-12) {
          const double factor = 2.0 * Sfm / denom *
            (((double)F[i] - mF) - (Sfm / Smm) * (Mw[i] - mM));
          fx[i] += wk * factor * gFx[i];
          fy[i] += wk * factor * gFy[i];
          fz[i] += wk * factor * gFz[i];
          accCC += (Sfm * Sfm) / denom; ++cnt;
        }
      }
      blockAcc[b] = accCC; blockCnt[b] = cnt;
    }
  }
};

// Mean-squares force, accumulated (weighted) into the shared force field; the
// SSE is reduced per block of kCostBlock voxels for the cost trace.
struct MSForceWorker : public TinyParallel::Worker {
  const reg_real* F = nullptr;
  const double *Mw = nullptr, *gFx = nullptr, *gFy = nullptr, *gFz = nullptr;
  reg_real *fx = nullptr, *fy = nullptr, *fz = nullptr;
  double wk = 0.0;
  const reg_real* mask = nullptr;   // optional fixed mask: skip force where 0
  std::size_t N = 0;
  double* blockSse = nullptr;       // per block: sum of the squared differences

  void operator()(std::size_t bBegin, std::size_t bEnd) override {
    for (std::size_t b = bBegin; b < bEnd; ++b) {
      const std::size_t begin = b * kCostBlock;
      const std::size_t end = std::min(N, begin + kCostBlock);
      double sse = 0.0;
      for (std::size_t i = begin; i < end; ++i) {
        if (mask && mask[i] == reg_real(0)) continue;   // outside fixed mask
        const double d = (double)F[i] - Mw[i];
        fx[i] += wk * d * gFx[i];
        fy[i] += wk * d * gFy[i];
        fz[i] += wk * d * gFz[i];
        sse += d * d;
      }
      blockSse[b] = sse;
    }
  }
};

// Maximum field magnitude (squared) over a reg_real field, for the step cap and
// the verbose field_max readout.
struct MaxNormReducer : public TinyParallel::Worker {
  const reg_real *fx, *fy, *fz;
  double maxn2;

  MaxNormReducer(const reg_real* fx_, const reg_real* fy_, const reg_real* fz_)
    : fx(fx_), fy(fy_), fz(fz_), maxn2(0.0) {}
  MaxNormReducer(MaxNormReducer& o, TinyParallel::Split)
    : fx(o.fx), fy(o.fy), fz(o.fz), maxn2(0.0) {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      // accumulate in reg_real (as the serial code did, force being reg_real)
      // so the step-cap norm, hence the field update, is bit-identical
      const reg_real ax = fx[i], ay = fy[i], az = fz[i];
      const double n2 = ax * ax + ay * ay + az * az;
      if (n2 > maxn2) maxn2 = n2;
    }
  }
  void join(const MaxNormReducer& rhs) { if (rhs.maxn2 > maxn2) maxn2 = rhs.maxn2; }
};

// cur += scale * force.
struct FieldUpdateWorker : public TinyParallel::Worker {
  reg_real *cx, *cy, *cz;
  const reg_real *fx, *fy, *fz;
  double scale;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      const double sx = scale * fx[i], sy = scale * fy[i], sz = scale * fz[i];
      cx[i] += sx; cy[i] += sy; cz[i] += sz;
    }
  }
};

// ---------------------------------------------------------------------------
// Channel combination of the multi-metric of ITK (the one ANTs uses):
// derivative = sum_k w_k d_k / |d_k| * sum_k |d_k|, so only the weights set the
// share of each channel. Sums of squares are reduced over fixed blocks of
// kCostBlock voxels, combined in block order (identical for any thread count).
// ---------------------------------------------------------------------------
template <typename T>
struct SqNormBlockWorker : public TinyParallel::Worker {
  const T *x = nullptr, *y = nullptr, *z = nullptr;
  std::size_t N = 0, block = 0;
  double* out = nullptr;
  void operator()(std::size_t bBegin, std::size_t bEnd) override {
    for (std::size_t b = bBegin; b < bEnd; ++b) {
      const std::size_t begin = b * block, end = std::min(N, begin + block);
      double s = 0.0;
      for (std::size_t i = begin; i < end; ++i) {
        const double a = x[i], c = y[i], d = z[i];
        s += a * a + c * c + d * d;
      }
      out[b] = s;
    }
  }
};

template <typename T>
double blockSumSquares(const T* x, const T* y, const T* z, std::size_t N,
                       std::vector<double>& blocks) {
  SqNormBlockWorker<T> nw;
  nw.x = x; nw.y = y; nw.z = z; nw.N = N; nw.block = kCostBlock; nw.out = blocks.data();
  TinyParallel::parallelFor(0, blocks.size(), nw, 1);
  double s = 0.0;
  for (std::size_t b = 0; b < blocks.size(); ++b) s += blocks[b];
  return s;
}

// dst += a * src (per voxel, independent)
struct AxpyFieldWorker : public TinyParallel::Worker {
  reg_real *dx, *dy, *dz;
  const reg_real *sx, *sy, *sz;
  double a;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      dx[i] += (reg_real)(a * sx[i]); dy[i] += (reg_real)(a * sy[i]); dz[i] += (reg_real)(a * sz[i]);
    }
  }
};

// f *= a (per voxel, independent)
struct ScaleFieldWorker : public TinyParallel::Worker {
  reg_real *fx, *fy, *fz;
  double a;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      fx[i] = (reg_real)(a * fx[i]); fy[i] = (reg_real)(a * fy[i]); fz[i] = (reg_real)(a * fz[i]);
    }
  }
};

// ---------------------------------------------------------------------------
// Mattes mutual information for the deformable stage (the MI metric of ANTs,
// ITK's MattesMutualInformationImageToImageMetricv4 with a dense transform).
// The joint histogram has `bins` bins per axis, 2 of them padding at each end;
// fixed intensities fall in one bin (a zero-order window), moving intensities
// are spread by a cubic B-spline Parzen window. The force at a sampled voxel is
//   -(1 / Δm) sum_κ β3'(κ - m̃) log( p(ι, κ) / p_F(ι) ) * grad M,
// ITK's pRatio, without its 1 / N (a constant that the normalized channel
// combination, or with one channel the step cap, removes). The histogram and
// the forces run serially over the samples, so the result does not depend on
// the thread count.
// ---------------------------------------------------------------------------
const int kMIPadding = 2;

inline double bspline3(double t) {
  const double a = std::abs(t);
  if (a < 1.0) return (4.0 - 6.0 * a * a + 3.0 * a * a * a) / 6.0;
  if (a < 2.0) { const double b = 2.0 - a; return b * b * b / 6.0; }
  return 0.0;
}

inline double bspline3Deriv(double t) {
  const double a = std::abs(t);
  double d;
  if (a < 1.0) {
    d = -2.0 * a + 1.5 * a * a;
  } else if (a < 2.0) {
    const double b = 2.0 - a; d = -0.5 * b * b;
  } else {
    return 0.0;
  }
  return (t < 0.0) ? -d : d;
}

struct MIRange {
  double fmin = 0.0, fbin = 1.0, mmin = 0.0, mbin = 1.0;
  int bins = 32;

  void set(double f0, double f1, double m0, double m1, int nbins) {
    bins = nbins;
    const double inner = (double)(bins - 2 * kMIPadding);
    fmin = f0; fbin = (f1 > f0) ? (f1 - f0) / inner : 1.0;
    mmin = m0; mbin = (m1 > m0) ? (m1 - m0) / inner : 1.0;
  }
  inline int fixedBin(double f) const {
    int b = (int)std::floor((f - fmin) / fbin) + kMIPadding;
    if (b < kMIPadding) b = kMIPadding;
    if (b > bins - kMIPadding - 1) b = bins - kMIPadding - 1;
    return b;
  }
  // continuous moving position m̃ and the first of its 4 Parzen bins (c - 1)
  inline void movingPos(double m, double& mt, int& c) const {
    mt = (m - mmin) / mbin + kMIPadding;
    if (mt < kMIPadding) mt = kMIPadding;
    if (mt > bins - kMIPadding) mt = bins - kMIPadding;
    c = (int)std::floor(mt);
    if (c < kMIPadding) c = kMIPadding;
    if (c > bins - kMIPadding - 1) c = bins - kMIPadding - 1;
  }
};

// The mutual information of the valid samples, and its force accumulated as
// fx[s] += wk * dMI/dm * grad M at those samples. `joint` and `logRatio` are
// caller-owned scratch (bins * bins).
double mattesForce(const reg_real* F, const double* Mw, const uint8_t* valid,
                   const double* gx, const double* gy, const double* gz,
                   const std::vector<std::size_t>& samples, const MIRange& R,
                   double wk, reg_real* fx, reg_real* fy, reg_real* fz,
                   std::vector<double>& joint, std::vector<double>& logRatio) {
  const int B = R.bins;
  joint.assign((std::size_t)B * B, 0.0);
  double total = 0.0;
  for (std::size_t s : samples) {
    if (!valid[s]) continue;
    const int iota = R.fixedBin(F[s]);
    double mt; int c; R.movingPos(Mw[s], mt, c);
    for (int kap = c - 1; kap <= c + 2; ++kap) {
      const double b = bspline3((double)kap - mt);
      joint[(std::size_t)iota * B + kap] += b;
      total += b;
    }
  }
  if (!(total > 0.0)) return 0.0;
  std::vector<double> pF(B, 0.0), pM(B, 0.0);
  for (int i = 0; i < B; ++i) {
    for (int k = 0; k < B; ++k) {
      double& p = joint[(std::size_t)i * B + k];
      p /= total;
      pF[i] += p; pM[k] += p;
    }
  }
  // the ratio is taken to the fixed marginal, unlike ITK's ComputeResults
  // (log(p / p_M), the exact derivative of MI): the fixed term cancels over
  // the Parzen window, so the force is the gradient of the negative joint
  // entropy. The moving-entropy term that this leaves out kept pushing aligned
  // images apart in the one-sided update (2.3 voxels on a phantom where
  // 'ANTs' moved 0.004, and 0.3 with this ratio)
  double mi = 0.0;
  logRatio.assign((std::size_t)B * B, 0.0);
  for (int i = 0; i < B; ++i) {
    for (int k = 0; k < B; ++k) {
      const double p = joint[(std::size_t)i * B + k];
      if (p > 1e-16) {
        mi += p * std::log(p / (pF[i] * pM[k]));
        logRatio[(std::size_t)i * B + k] = std::log(p / pF[i]);
      }
    }
  }
  for (std::size_t s : samples) {
    if (!valid[s]) continue;
    const int iota = R.fixedBin(F[s]);
    double mt; int c; R.movingPos(Mw[s], mt, c);
    double d = 0.0;
    for (int kap = c - 1; kap <= c + 2; ++kap) {
      d += bspline3Deriv((double)kap - mt) * logRatio[(std::size_t)iota * B + kap];
    }
    const double factor = -d / R.mbin;
    fx[s] += wk * factor * gx[s];
    fy[s] += wk * factor * gy[s];
    fz[s] += wk * factor * gz[s];
  }
  return mi;
}

// Sample a moving image at a RAS point with the requested output interpolation:
// 0 = nearest (label-preserving), 1 = trilinear, 2 = cubic B-spline. Used only
// for the final warped output, the optimizer's metric sampling stays trilinear
// so every channel produces smooth gradients regardless of this choice.
inline double sampleWithCode(const RegImage3& img, const Vector3d& qras, int code) {
  Vector4d ph(qras[0], qras[1], qras[2], 1.0);
  Vector4d v = img.ras2vox * ph;
  const double cx = v[0], cy = v[1], cz = v[2];
  double val = 0.0;
  if (code == 0) {                 // nearest neighbor
    const long ix = std::lround(cx), iy = std::lround(cy), iz = std::lround(cz);
    if (ix < 0 || iy < 0 || iz < 0 || ix >= img.nx || iy >= img.ny || iz >= img.nz)
      return 0.0;
    return static_cast<double>(
        img.data[ix + static_cast<std::size_t>(img.nx) * (iy + static_cast<std::size_t>(img.ny) * iz)]);
  } else if (code == 2) {          // cubic B-spline
    bsplineSample(img.data, img.nx, img.ny, img.nz, cx, cy, cz, val);
  } else {                         // trilinear
    trilinearSample(img.data, img.nx, img.ny, img.nz, cx, cy, cz, val);
  }
  return val;
}

// Final per-channel warp of every moving channel onto the full fixed grid,
// through aff o (id + full_disp), each with its own output interpolation.
// Parallelized over the fixed grid; writes to the raw data of pre-allocated
// output vectors at disjoint indices (no R API touched in the worker).
struct OutputWarpWorker : public TinyParallel::Worker {
  Matrix4d fV2R, aff;
  const reg_real *dispx, *dispy, *dispz;
  int fnx; std::size_t szf;
  const RegImage3* movFulls;
  const int* codes;
  int K;
  double* const* outs;             // K output buffers

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      int i, j, k; decodeIndex(o, fnx, szf, i, j, k);
      Vector4d ph(i, j, k, 1.0);
      Vector3d x = (fV2R * ph).head<3>();
      Vector4d xdh(x[0] + dispx[o], x[1] + dispy[o], x[2] + dispz[o], 1.0);
      Vector3d q = (aff * xdh).head<3>();
      for (int kk = 0; kk < K; ++kk)
        outs[kk][o] = sampleWithCode(movFulls[kk], q, codes[kk]);
    }
  }
};

// Trilinear "splat" (scatter): the adjoint of trilinearSample's gather. Adds
// w * vec into the field at continuous level-voxel coords (cx,cy,cz), spread to
// the 8 surrounding voxels with the same clamped indices/weights the gather
// uses. Out-of-bounds coords contribute nothing (matching a false gather). Used
// by the cortical-landmark force, which is sparse, so this runs serially.
inline void splatTrilinear(Field& f, double cx, double cy, double cz,
                           double w, const Vector3d& vec) {
  if (cx < 0.0 || cy < 0.0 || cz < 0.0 ||
      cx > (double)(f.nx - 1) || cy > (double)(f.ny - 1) || cz > (double)(f.nz - 1))
    return;
  const int x0 = (int)std::floor(cx), y0 = (int)std::floor(cy), z0 = (int)std::floor(cz);
  const int x1 = (x0 < f.nx - 1) ? x0 + 1 : x0;
  const int y1 = (y0 < f.ny - 1) ? y0 + 1 : y0;
  const int z1 = (z0 < f.nz - 1) ? z0 + 1 : z0;
  const double fx = cx - x0, fy = cy - y0, fz = cz - z0;
  const std::size_t sy = (std::size_t)f.nx, sz = (std::size_t)f.nx * f.ny;
  const int xs[2] = {x0, x1}, ys[2] = {y0, y1}, zs[2] = {z0, z1};
  const double wxw[2] = {1.0 - fx, fx}, wyw[2] = {1.0 - fy, fy}, wzw[2] = {1.0 - fz, fz};
  for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b)
      for (int c = 0; c < 2; ++c) {
        const double cw = w * wxw[a] * wyw[b] * wzw[c];
        if (cw == 0.0) continue;     // skips the duplicated clamped corner
        const std::size_t o = (std::size_t)xs[a] + sy * ys[b] + sz * zs[c];
        f.x[o] += (reg_real)(cw * vec[0]);
        f.y[o] += (reg_real)(cw * vec[1]);
        f.z[o] += (reg_real)(cw * vec[2]);
      }
}

// The convergence rule of the deformable stages of ANTs (ITK's
// WindowConvergenceMonitoringFunction, as SyNImageRegistrationMethod uses it).
// Every energy E (the cost of the first channel) is added to a running sum
// T = sum |E| over the level. Once W energies exist, the last W, divided by T
// and placed at u = 0, 1/(W-1), ..., 1 (the last one clamped to 1 - 1e-4, as
// in ITK), are fitted by ITK's one-level, order-1 B-spline approximation (two
// control points phi_0 and phi_1), and the convergence value is the fitted drop
// across the window, phi_0 - phi_1 = sum_n c_n E_n / T. It is +Inf until the
// window is full (no stop before W iterations), and 0 when the quotient is not
// finite (ITK replaces such control points by 0). The weights c_n are computed
// when the window first fills, so a window longer than the level costs nothing.
struct ConvergenceMonitor {
  int W;
  std::vector<double> c;            // fit weights
  std::deque<double> window;        // the last W energies
  double total = 0.0;               // sum of |E| since the level started

  explicit ConvergenceMonitor(int W_) : W(W_) {}

  void computeWeights() {
    std::vector<double> u(W);
    for (int n = 0; n < W; ++n) u[n] = (double)n / (double)(W - 1);
    u[W - 1] = 1.0 - 1e-4;
    double A = 0.0, B = 0.0;
    for (int n = 0; n < W; ++n) {
      A += (1.0 - u[n]) * (1.0 - u[n]);
      B += u[n] * u[n];
    }
    c.assign(W, 0.0);
    for (int n = 0; n < W; ++n) {
      const double a = 1.0 - u[n], b = u[n], S = a * a + b * b;
      c[n] = a * a * a / (S * A) - b * b * b / (S * B);
    }
  }

  double add(double e) {
    window.push_back(e);
    if (window.size() > (std::size_t)W) window.pop_front();
    total += std::abs(e);
    if (window.size() < (std::size_t)W) return R_PosInf;
    if (c.empty()) computeWeights();
    double s = 0.0;
    for (int n = 0; n < W; ++n) s += c[n] * window[n];
    const double v = s / total;
    return std::isfinite(v) ? v : 0.0;
  }
};

// ===========================================================================
// Symmetric SyN (as SyNImageRegistrationMethod of ITK, used by ANTs). Two
// displacement fields live on the fixed level grid, the virtual ("middle")
// domain: uF maps the middle to the fixed image and uM the middle to the moving
// image before the affine. Every iteration warps both images into the middle,
// computes the force of each side from the gradient of its own warped image
// ((I + grad u)^T grad I), smooths and caps each update, and composes it with
// its field (u <- delta + u(x + delta)), which keeps the maps one-to-one. The
// result is A(phi_M(phi_F^-1(x))): the inverse of uF composed with uM.
// ===========================================================================

// The fixed channel (smoothed, full grid) sampled at x + uF(x), with its RAS
// gradient; `valid` (optional) marks samples inside the image and the gate.
struct FixedWarpWorker : public TinyParallel::Worker {
  const RegImage3* img = nullptr;
  const double *bx = nullptr, *by = nullptr, *bz = nullptr;
  const reg_real *ux = nullptr, *uy = nullptr, *uz = nullptr;
  double* vD = nullptr; reg_real* vR = nullptr;
  double *gx = nullptr, *gy = nullptr, *gz = nullptr;
  uint8_t* valid = nullptr;
  const uint8_t* gate = nullptr;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      const Vector3d p(bx[o] + ux[o], by[o] + uy[o], bz[o] + uz[o]);
      double v; Vector3d g;
      const bool inside = img->sampleRas(p, v, &g);
      vD[o] = v; vR[o] = (reg_real)v;
      gx[o] = g[0]; gy[o] = g[1]; gz[o] = g[2];
      if (valid) valid[o] = (inside && (!gate || gate[o])) ? 1 : 0;
    }
  }
};

// The fixed mask sampled at x + uF(x): the gate of the forces in the middle.
struct MaskGateWorker : public TinyParallel::Worker {
  const RegImage3* mask = nullptr;
  const double *bx = nullptr, *by = nullptr, *bz = nullptr;
  const reg_real *ux = nullptr, *uy = nullptr, *uz = nullptr;
  uint8_t* gate = nullptr;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      const Vector3d p(bx[o] + ux[o], by[o] + uy[o], bz[o] + uz[o]);
      double mv; mask->sampleRas(p, mv);
      gate[o] = (mv > 0.5) ? 1 : 0;
    }
  }
};

// g <- (I + grad u)^T g, the gradient of the warped image in the middle, with
// grad u from central differences of the field (one-sided at the borders).
struct JacobianTWorker : public TinyParallel::Worker {
  const reg_real *ux = nullptr, *uy = nullptr, *uz = nullptr;
  int nx = 0, ny = 0, nz = 0;
  Matrix3d ras2voxR = Matrix3d::Identity();   // level grid
  double *gx = nullptr, *gy = nullptr, *gz = nullptr;
  void operator()(std::size_t begin, std::size_t end) override {
    const std::size_t sy = (std::size_t)nx, sz = (std::size_t)nx * ny;
    const int n[3] = {nx, ny, nz};
    const std::size_t st[3] = {1, sy, sz};
    const reg_real* u[3] = {ux, uy, uz};
    for (std::size_t o = begin; o < end; ++o) {
      int i, j, k; decodeIndex(o, nx, sz, i, j, k);
      const int idx[3] = {i, j, k};
      double D[3][3];                       // d u_c / d index_a
      for (int a = 0; a < 3; ++a) {
        const int ia = idx[a];
        if (n[a] < 2) { for (int c = 0; c < 3; ++c) D[c][a] = 0.0; continue; }
        const std::size_t lo = (ia > 0) ? o - st[a] : o;
        const std::size_t hi = (ia < n[a] - 1) ? o + st[a] : o;
        const double den = (ia > 0 && ia < n[a] - 1) ? 2.0 : 1.0;
        for (int c = 0; c < 3; ++c) D[c][a] = ((double)u[c][hi] - (double)u[c][lo]) / den;
      }
      const double gc[3] = {gx[o], gy[o], gz[o]};
      double out[3];
      for (int jj = 0; jj < 3; ++jj) {
        // (grad u)_{c, jj} = sum_a D[c][a] * ras2vox(a, jj)
        double add = 0.0;
        for (int c = 0; c < 3; ++c) {
          double gu = 0.0;
          for (int a = 0; a < 3; ++a) gu += D[c][a] * ras2voxR(a, jj);
          add += gu * gc[c];
        }
        out[jj] = gc[jj] + add;
      }
      gx[o] = out[0]; gy[o] = out[1]; gz[o] = out[2];
    }
  }
};

// Local-CC forces of both sides in the middle (the derivative with respect to
// the fixed image swaps the roles of F and M), with the cost per block. Only
// the samples inside both images and both masks count (`valid`), as in ITK:
// a moving sample outside the moving mask is background, which must not pull
// the fixed half either.
struct SymCCWorker : public TinyParallel::Worker {
  const reg_real* F = nullptr; const double* M = nullptr;
  const double *sumF = nullptr, *sumFF = nullptr, *sumM = nullptr, *sumMM = nullptr, *sumFM = nullptr;
  const double *gFx = nullptr, *gFy = nullptr, *gFz = nullptr;
  const double *gMx = nullptr, *gMy = nullptr, *gMz = nullptr;
  reg_real *fFx = nullptr, *fFy = nullptr, *fFz = nullptr;
  reg_real *fMx = nullptr, *fMy = nullptr, *fMz = nullptr;
  double winN = 0.0, tau = 0.0, wk = 0.0;
  const uint8_t* valid = nullptr;
  std::size_t N = 0;
  double* blockAcc = nullptr; long* blockCnt = nullptr;
  void operator()(std::size_t bBegin, std::size_t bEnd) override {
    for (std::size_t b = bBegin; b < bEnd; ++b) {
      const std::size_t begin = b * kCostBlock, end = std::min(N, begin + kCostBlock);
      double acc = 0.0; long cnt = 0;
      for (std::size_t i = begin; i < end; ++i) {
        if (!valid[i]) continue;
        const double mF = sumF[i] / winN, mM = sumM[i] / winN;
        const double Sff = sumFF[i] - winN * mF * mF;
        const double Smm = sumMM[i] - winN * mM * mM;
        const double Sfm = sumFM[i] - winN * mF * mM;
        const double denom = Sff * Smm;
        if (Sff > tau && Smm > tau && denom > 1e-12) {
          const double dF = (double)F[i] - mF, dM = M[i] - mM;
          const double factorM = 2.0 * Sfm / denom * (dF - (Sfm / Smm) * dM);
          const double factorF = 2.0 * Sfm / denom * (dM - (Sfm / Sff) * dF);
          fMx[i] += wk * factorM * gMx[i]; fMy[i] += wk * factorM * gMy[i]; fMz[i] += wk * factorM * gMz[i];
          fFx[i] += wk * factorF * gFx[i]; fFy[i] += wk * factorF * gFy[i]; fFz[i] += wk * factorF * gFz[i];
          acc += (Sfm * Sfm) / denom; ++cnt;
        }
      }
      blockAcc[b] = acc; blockCnt[b] = cnt;
    }
  }
};

// Mean-squares forces of both sides in the middle (the samples of `valid`
// only, as SymCCWorker).
struct SymMSWorker : public TinyParallel::Worker {
  const reg_real* F = nullptr; const double* M = nullptr;
  const double *gFx = nullptr, *gFy = nullptr, *gFz = nullptr;
  const double *gMx = nullptr, *gMy = nullptr, *gMz = nullptr;
  reg_real *fFx = nullptr, *fFy = nullptr, *fFz = nullptr;
  reg_real *fMx = nullptr, *fMy = nullptr, *fMz = nullptr;
  double wk = 0.0;
  const uint8_t* valid = nullptr;
  std::size_t N = 0;
  double* blockSse = nullptr;
  void operator()(std::size_t bBegin, std::size_t bEnd) override {
    for (std::size_t b = bBegin; b < bEnd; ++b) {
      const std::size_t begin = b * kCostBlock, end = std::min(N, begin + kCostBlock);
      double sse = 0.0;
      for (std::size_t i = begin; i < end; ++i) {
        if (!valid[i]) continue;
        const double d = (double)F[i] - M[i];
        fMx[i] += wk * d * gMx[i]; fMy[i] += wk * d * gMy[i]; fMz[i] += wk * d * gMz[i];
        fFx[i] -= wk * d * gFx[i]; fFy[i] -= wk * d * gFy[i]; fFz[i] -= wk * d * gFz[i];
        sse += d * d;
      }
      blockSse[b] = sse;
    }
  }
};

// u_new(x) = delta(x) + u(x + delta(x)), u sampled trilinearly on its grid
// (coordinates clamped to the grid).
template <typename TU>
struct ComposeWorker : public TinyParallel::Worker {
  const double *bx = nullptr, *by = nullptr, *bz = nullptr;
  const reg_real *dx = nullptr, *dy = nullptr, *dz = nullptr;   // delta
  const TU *ux = nullptr, *uy = nullptr, *uz = nullptr;         // field
  Matrix4d ras2vox = Matrix4d::Identity();
  int nx = 0, ny = 0, nz = 0;
  reg_real *ox = nullptr, *oy = nullptr, *oz = nullptr;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      const double px = bx[o] + dx[o], py = by[o] + dy[o], pz = bz[o] + dz[o];
      const Vector4d v = ras2vox * Vector4d(px, py, pz, 1.0);
      const double cx = std::min(std::max(v[0], 0.0), (double)(nx - 1));
      const double cy = std::min(std::max(v[1], 0.0), (double)(ny - 1));
      const double cz = std::min(std::max(v[2], 0.0), (double)(nz - 1));
      double sx = 0.0, sy = 0.0, sz = 0.0;
      trilinearSample(ux, nx, ny, nz, cx, cy, cz, sx);
      trilinearSample(uy, nx, ny, nz, cx, cy, cz, sy);
      trilinearSample(uz, nx, ny, nz, cx, cy, cz, sz);
      ox[o] = (reg_real)(dx[o] + sx); oy[o] = (reg_real)(dy[o] + sy); oz[o] = (reg_real)(dz[o] + sz);
    }
  }
};

// Zero displacement on the outermost voxel layer, as the smoothing of the
// displacement fields in ITK's SyN ("make sure boundary does not move"): the
// border samples stay on the grid instead of falling outside it
void zeroFieldBoundary(Field& f) {
  const int nx = f.nx, ny = f.ny, nz = f.nz;
  const std::size_t sy = (std::size_t)nx, sz = (std::size_t)nx * ny;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      const bool edgeJK = (j == 0 || j == ny - 1 || k == 0 || k == nz - 1);
      for (int i = 0; i < nx; ++i) {
        if (edgeJK || i == 0 || i == nx - 1) {
          const std::size_t o = (std::size_t)i + sy * j + sz * k;
          f.x[o] = 0; f.y[o] = 0; f.z[o] = 0;
        }
      }
    }
  }
}

Rcpp::List registerSynSymmetric(
    const std::vector<Rcpp::NumericVector>& fixedVecs,
    const std::vector<Rcpp::NumericVector>& movingVecs,
    const int fnx, const int fny, const int fnz,
    const int mnx, const int mny, const int mnz,
    const Matrix4d& fV2R, const Matrix4d& mV2R, const Matrix4d& aff,
    const std::vector<int>& mcode, const std::vector<double>& w,
    const Rcpp::IntegerVector& shrinkFactors, const Rcpp::NumericVector& smoothingSigmas,
    const Rcpp::IntegerVector& iterations, const double gradStep,
    const Rcpp::NumericVector& flowSigma, const double totalSigma, const int ccRadius,
    const Rcpp::IntegerVector& interpCodes, const reg_real* fixedMaskData,
    const RegImage3* movMaskPtr, const std::vector<reg_real>& movMaskStore,
    const int miBins, const double miSampling, const bool normalizeChannels,
    const int seed, const bool monitorOn, const double convergenceThreshold,
    const int convergenceWindow, const bool verbose) {
  const Matrix3d affRT = aff.topLeftCorner<3, 3>().transpose();
  const int K = (int)fixedVecs.size();
  bool anyMI = false;
  for (int k = 0; k < K; ++k) if (mcode[k] == 2) anyMI = true;
  // the fixed mask on the full grid, sampled at x + uF(x) for the gate
  RegImage3 fixedMaskImg;
  if (fixedMaskData) fixedMaskImg = RegImage3(fixedMaskData, fnx, fny, fnz, fV2R);

  const int nLevels = shrinkFactors.size();
  std::vector<double> trace;
  Rcpp::IntegerVector iterationsRun(nLevels);
  Field dispF, dispM;
  bool haveField = false;

  for (int lev = 0; lev < nLevels; ++lev) {
    const int shrink = shrinkFactors[lev];
    const double sigma = (lev < smoothingSigmas.size()) ? smoothingSigmas[lev] : 0.0;
    const int maxIter = (lev < iterations.size()) ? iterations[lev] : 0;
    const double flowSigmaLev = (flowSigma.size() == 0) ? 0.0
        : (lev < flowSigma.size() ? flowSigma[lev] : flowSigma[flowSigma.size() - 1]);

    // smoothed channels on their full grids: both sides are sampled at warped
    // points; the fixed level images give the grid and the per-level statistics
    std::vector<std::vector<reg_real> > sFixedStore(K), sMovingStore(K);
    std::vector<FixedLevel> flvls; flvls.reserve(K);
    std::vector<RegImage3> fixs, movs; fixs.reserve(K); movs.reserve(K);
    for (int k = 0; k < K; ++k) {
      sFixedStore[k]  = gaussianSmooth(&fixedVecs[k][0], fnx, fny, fnz, sigma);
      sMovingStore[k] = gaussianSmooth(&movingVecs[k][0], mnx, mny, mnz, sigma);
      flvls.push_back(buildFixedLevel(sFixedStore[k].data(), fnx, fny, fnz, fV2R, shrink));
      fixs.push_back(RegImage3(sFixedStore[k].data(), fnx, fny, fnz, fV2R));
      movs.push_back(RegImage3(sMovingStore[k].data(), mnx, mny, mnz, mV2R));
    }
    const FixedLevel& flvl0 = flvls[0];
    const int lnx = flvl0.nx, lny = flvl0.ny, lnz = flvl0.nz;
    const std::size_t N = flvl0.size();
    const std::size_t szl = (std::size_t)lnx * (std::size_t)lny;
    const Matrix4d levelV2R = flvl0.vox2ras;
    const Matrix4d levelR2V = levelV2R.inverse();
    const Matrix3d levelR2VR = levelR2V.topLeftCorner<3, 3>();

    Field curF, curM;
    curF.nx = lnx; curF.ny = lny; curF.nz = lnz; curF.vox2ras = levelV2R; curF.alloc(lnx, lny, lnz);
    curM.nx = lnx; curM.ny = lny; curM.nz = lnz; curM.vox2ras = levelV2R; curM.alloc(lnx, lny, lnz);
    if (haveField) {
      resampleField(dispF, curF); resampleField(dispM, curM);
      zeroFieldBoundary(curF); zeroFieldBoundary(curM);
    }

    const double levelSpacing = std::cbrt(std::abs(levelV2R.topLeftCorner<3, 3>().determinant()));

    if (verbose) {
      Rcpp::Rcout << "[SyN, symmetric] level " << (lev + 1) << "/" << nLevels
                  << " (shrink=" << shrink << ", sigma=" << std::fixed << std::setprecision(1) << sigma
                  << ", flow=" << std::setprecision(1) << flowSigmaLev << ", channels=" << K
                  << "): max " << maxIter << " iterations" << std::endl;
    }

    if (maxIter > 0) {
      const double winN = std::pow(2.0 * ccRadius + 1.0, 3);
      std::vector<double> bx(N), by(N), bz(N);
      {
        BaseGridWorker bw;
        bw.v2r = levelV2R; bw.nx = lnx; bw.sz = szl;
        bw.bx = bx.data(); bw.by = by.data(); bw.bz = bz.data();
        TinyParallel::parallelFor(0, N, bw, 1024);
      }
      // the static level mask (for the MI samples and the per-level statistics)
      std::vector<reg_real> maskLev;
      if (fixedMaskData) {
        FixedLevel ml = buildFixedLevel(fixedMaskData, fnx, fny, fnz, fV2R, shrink);
        maskLev.resize(N);
        for (std::size_t i = 0; i < N; ++i) maskLev[i] = (ml.data[i] > reg_real(0.5)) ? reg_real(1) : reg_real(0);
      }
      const reg_real* maskLevPtr = maskLev.empty() ? nullptr : maskLev.data();

      // per channel: the CC threshold, the MI ranges, the noise-floor scale
      std::vector<double> ccTau(K, 0.0), chScale(K, 0.0);
      std::vector<MIRange> miRange(K), miRangeSwap(K);
      for (int k = 0; k < K; ++k) {
        const std::vector<reg_real>& Fk = flvls[k].data;
        double ga = 0.0, gv = 0.0;
        for (std::size_t i = 0; i < N; ++i) ga += Fk[i];
        ga /= (double)N;
        for (std::size_t i = 0; i < N; ++i) { const double d = Fk[i] - ga; gv += d * d; }
        gv /= (double)N;
        ccTau[k] = 1e-3 * gv * winN;
        double s1 = 0.0, s2 = 0.0, f0 = R_PosInf, f1 = R_NegInf; std::size_t n = 0;
        for (std::size_t i = 0; i < N; ++i) {
          if (maskLevPtr && maskLevPtr[i] == reg_real(0)) continue;
          const double v = Fk[i]; s1 += v; s2 += v * v; ++n;
          if (v < f0) f0 = v;
          if (v > f1) f1 = v;
        }
        double sd = 0.0;
        if (n > 1) { const double m = s1 / (double)n; sd = std::sqrt(std::max(0.0, s2 / (double)n - m * m)); }
        if (mcode[k] == 2) {
          const std::vector<reg_real>& Mk = sMovingStore[k];
          double m0 = R_PosInf, m1 = R_NegInf;
          for (std::size_t i = 0; i < Mk.size(); ++i) {
            if (movMaskPtr && movMaskStore[i] <= reg_real(0.5)) continue;
            const double v = Mk[i];
            if (v < m0) m0 = v;
            if (v > m1) m1 = v;
          }
          if (!std::isfinite(f0)) { f0 = 0.0; f1 = 1.0; }
          if (!std::isfinite(m0)) { m0 = 0.0; m1 = 1.0; }
          miRange[k].set(f0, f1, m0, m1, miBins);       // moving side: cubic on M
          miRangeSwap[k].set(m0, m1, f0, f1, miBins);   // fixed side: cubic on F
          chScale[k] = 1.0 / std::min(miRange[k].mbin, miRange[k].fbin);
        } else {
          chScale[k] = (mcode[k] == 0) ? sd : ((sd > 0.0) ? 2.0 / (sd * winN) : 0.0);
        }
      }
      std::vector<std::size_t> miSamples;
      if (anyMI) {
        std::mt19937 rng((unsigned int)seed + 7919u * (unsigned int)(lev + 1));
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        for (std::size_t i = 0; i < N; ++i) {
          if (maskLevPtr && maskLevPtr[i] == reg_real(0)) continue;
          if (miSampling >= 1.0 || unif(rng) < miSampling) miSamples.push_back(i);
        }
      }

      // scratch
      std::vector<double> FmD(N), Mw(N), gFx(N), gFy(N), gFz(N), gMx(N), gMy(N), gMz(N);
      std::vector<reg_real> FmR(N), MwR, gscratch(N);
      std::vector<double> sumF(N), sumFF(N), sumM(N), sumMM(N), sumFM(N);
      std::vector<double> scrM(N), scrMM(N), scrFM(N), FF(N), scrF(N);
      // validF: the fixed sample is inside the fixed image (and mask);
      // validM: the moving sample is inside the moving image and mask; the
      // forces and costs of both halves take the samples valid on both sides
      std::vector<uint8_t> validF(N, 0), validM(N, 0), validB(N, 0), gate;
      if (anyMI) MwR.assign(N, 0);
      if (fixedMaskData) gate.assign(N, 0);
      Field forceF, forceM, chF, chM, tmp;
      Field* fl[5] = {&forceF, &forceM, &chF, &chM, &tmp};
      for (int t = 0; t < 5; ++t) {
        if (t >= 2 && t <= 3 && !normalizeChannels) continue;
        fl[t]->nx = lnx; fl[t]->ny = lny; fl[t]->nz = lnz; fl[t]->vox2ras = levelV2R; fl[t]->alloc(lnx, lny, lnz);
      }
      const std::size_t nCostBlocks = (N + kCostBlock - 1) / kCostBlock;
      std::vector<double> blockAcc(nCostBlocks), blockSse(nCostBlocks), normBlocks(nCostBlocks);
      std::vector<long> blockCnt(nCostBlocks);
      std::vector<double> miJoint, miLog;

      ConvergenceMonitor monitor(monitorOn ? convergenceWindow : 0);
      int itRun = 0, stopCode = 0;
      double lastCost = 0.0, lastConv = R_PosInf;

      for (int it = 0; it < maxIter; ++it) {
        forceF.zero(); forceM.zero();
        double metricVal = 0.0, energy0 = 0.0, magF = 0.0, magM = 0.0;

        if (fixedMaskData) {
          MaskGateWorker mg;
          mg.mask = &fixedMaskImg; mg.bx = bx.data(); mg.by = by.data(); mg.bz = bz.data();
          mg.ux = curF.x.data(); mg.uy = curF.y.data(); mg.uz = curF.z.data();
          mg.gate = gate.data();
          TinyParallel::parallelFor(0, N, mg, 256);
        }

        for (int k = 0; k < K; ++k) {
          const double wk = w[k];
          Field& tF = normalizeChannels ? chF : forceF;
          Field& tM = normalizeChannels ? chM : forceM;
          const double wForce = normalizeChannels ? 1.0 : wk;
          if (normalizeChannels) { chF.zero(); chM.zero(); }

          // both sides into the middle
          {
            FixedWarpWorker fw;
            fw.img = &fixs[k]; fw.bx = bx.data(); fw.by = by.data(); fw.bz = bz.data();
            fw.ux = curF.x.data(); fw.uy = curF.y.data(); fw.uz = curF.z.data();
            fw.vD = FmD.data(); fw.vR = FmR.data(); fw.gx = gFx.data(); fw.gy = gFy.data(); fw.gz = gFz.data();
            fw.valid = validF.data();
            fw.gate = fixedMaskData ? gate.data() : nullptr;
            TinyParallel::parallelFor(0, N, fw, 256);
          }
          {
            WarpWorker ww;
            ww.mov = &movs[k]; ww.aff = aff; ww.affRT = affRT;
            ww.movMask = movMaskPtr; ww.warpMask = nullptr;
            ww.bx = bx.data(); ww.by = by.data(); ww.bz = bz.data();
            ww.dispx = curM.x.data(); ww.dispy = curM.y.data(); ww.dispz = curM.z.data();
            ww.Mw = Mw.data(); ww.gFx = gMx.data(); ww.gFy = gMy.data(); ww.gFz = gMz.data();
            ww.valid = validM.data();
            TinyParallel::parallelFor(0, N, ww, 256);
          }
          for (std::size_t i = 0; i < N; ++i) validB[i] = (validF[i] && validM[i]) ? 1 : 0;
          // gradients of the warped images in the middle
          {
            JacobianTWorker jw;
            jw.nx = lnx; jw.ny = lny; jw.nz = lnz; jw.ras2voxR = levelR2VR;
            jw.ux = curF.x.data(); jw.uy = curF.y.data(); jw.uz = curF.z.data();
            jw.gx = gFx.data(); jw.gy = gFy.data(); jw.gz = gFz.data();
            TinyParallel::parallelFor(0, N, jw, 1024);
            jw.ux = curM.x.data(); jw.uy = curM.y.data(); jw.uz = curM.z.data();
            jw.gx = gMx.data(); jw.gy = gMy.data(); jw.gz = gMz.data();
            TinyParallel::parallelFor(0, N, jw, 1024);
          }

          if (mcode[k] == 1) {
            boxSumInto(FmR.data(), lnx, lny, lnz, ccRadius, sumF, scrF);
            for (std::size_t i = 0; i < N; ++i) FF[i] = FmD[i] * FmD[i];
            boxSumInto(FF.data(), lnx, lny, lnz, ccRadius, sumFF, scrF);
            fusedLocalCCSums(Mw.data(), FmR.data(), lnx, lny, lnz, ccRadius,
                             sumM, sumMM, sumFM, scrM, scrMM, scrFM);
            SymCCWorker cc;
            cc.F = FmR.data(); cc.M = Mw.data();
            cc.sumF = sumF.data(); cc.sumFF = sumFF.data(); cc.sumM = sumM.data();
            cc.sumMM = sumMM.data(); cc.sumFM = sumFM.data();
            cc.gFx = gFx.data(); cc.gFy = gFy.data(); cc.gFz = gFz.data();
            cc.gMx = gMx.data(); cc.gMy = gMy.data(); cc.gMz = gMz.data();
            cc.fFx = tF.x.data(); cc.fFy = tF.y.data(); cc.fFz = tF.z.data();
            cc.fMx = tM.x.data(); cc.fMy = tM.y.data(); cc.fMz = tM.z.data();
            cc.winN = winN; cc.tau = ccTau[k]; cc.wk = wForce; cc.valid = validB.data(); cc.N = N;
            cc.blockAcc = blockAcc.data(); cc.blockCnt = blockCnt.data();
            TinyParallel::parallelFor(0, nCostBlocks, cc, 1);
            double acc = 0.0; long cnt = 0;
            for (std::size_t b = 0; b < nCostBlocks; ++b) { acc += blockAcc[b]; cnt += blockCnt[b]; }
            const double ek = (cnt > 0) ? -(acc / cnt) : 0.0;
            if (k == 0) energy0 = ek;
            metricVal += wk * ek;
          } else if (mcode[k] == 2) {
            for (std::size_t i = 0; i < N; ++i) MwR[i] = (reg_real)Mw[i];
            const double mi = mattesForce(FmR.data(), Mw.data(), validB.data(),
                                          gMx.data(), gMy.data(), gMz.data(), miSamples,
                                          miRange[k], wForce, tM.x.data(), tM.y.data(), tM.z.data(),
                                          miJoint, miLog);
            mattesForce(MwR.data(), FmD.data(), validB.data(), gFx.data(), gFy.data(), gFz.data(),
                        miSamples, miRangeSwap[k], wForce, tF.x.data(), tF.y.data(), tF.z.data(),
                        miJoint, miLog);
            const double ek = -mi;
            if (k == 0) energy0 = ek;
            metricVal += wk * ek;
          } else {
            SymMSWorker ms;
            ms.F = FmR.data(); ms.M = Mw.data();
            ms.gFx = gFx.data(); ms.gFy = gFy.data(); ms.gFz = gFz.data();
            ms.gMx = gMx.data(); ms.gMy = gMy.data(); ms.gMz = gMz.data();
            ms.fFx = tF.x.data(); ms.fFy = tF.y.data(); ms.fFz = tF.z.data();
            ms.fMx = tM.x.data(); ms.fMy = tM.y.data(); ms.fMz = tM.z.data();
            ms.wk = wForce; ms.valid = validB.data(); ms.N = N; ms.blockSse = blockSse.data();
            TinyParallel::parallelFor(0, nCostBlocks, ms, 1);
            double sse = 0.0;
            for (std::size_t b = 0; b < nCostBlocks; ++b) sse += blockSse[b];
            const double ek = sse / N;
            if (k == 0) energy0 = ek;
            metricVal += wk * ek;
          }

          if (normalizeChannels) {
            // each side: add w_k d_k / |d_k| above the noise floor
            Field* chs[2] = {&chF, &chM};
            Field* dst[2] = {&forceF, &forceM};
            const double* gs[2][3] = {{gFx.data(), gFy.data(), gFz.data()}, {gMx.data(), gMy.data(), gMz.data()}};
            double* mags[2] = {&magF, &magM};
            for (int side = 0; side < 2; ++side) {
              const double nk = std::sqrt(blockSumSquares(chs[side]->x.data(), chs[side]->y.data(),
                                                          chs[side]->z.data(), N, normBlocks));
              const double gk = std::sqrt(blockSumSquares(gs[side][0], gs[side][1], gs[side][2], N, normBlocks));
              if (nk > 1e-30 && nk > 1e-6 * gk * chScale[k]) {
                *mags[side] += nk;
                AxpyFieldWorker aw;
                aw.dx = dst[side]->x.data(); aw.dy = dst[side]->y.data(); aw.dz = dst[side]->z.data();
                aw.sx = chs[side]->x.data(); aw.sy = chs[side]->y.data(); aw.sz = chs[side]->z.data();
                aw.a = wk / nk;
                TinyParallel::parallelFor(0, N, aw, 4096);
              }
            }
          }
        } // channels

        // each side: (rescale,) smooth, cap at grad_step * spacing, compose
        Field* forces[2] = {&forceF, &forceM};
        Field* curs[2] = {&curF, &curM};
        const double mags2[2] = {magF, magM};
        double maxnSide[2] = {0.0, 0.0};
        for (int side = 0; side < 2; ++side) {
          Field& f = *forces[side];
          if (normalizeChannels && mags2[side] > 0.0) {
            ScaleFieldWorker sw;
            sw.fx = f.x.data(); sw.fy = f.y.data(); sw.fz = f.z.data(); sw.a = mags2[side];
            TinyParallel::parallelFor(0, N, sw, 4096);
          }
          smoothFieldInto(f, flowSigmaLev, gscratch.data());
          zeroFieldBoundary(f);
          MaxNormReducer mn(f.x.data(), f.y.data(), f.z.data());
          TinyParallel::parallelReduce(0, N, mn, 4096);
          const double maxn = std::sqrt(mn.maxn2);
          maxnSide[side] = maxn;
          if (maxn <= 1e-12) continue;
          {
            ScaleFieldWorker sw;
            sw.fx = f.x.data(); sw.fy = f.y.data(); sw.fz = f.z.data(); sw.a = gradStep * levelSpacing / maxn;
            TinyParallel::parallelFor(0, N, sw, 4096);
          }
          Field& u = *curs[side];
          {
            ComposeWorker<reg_real> cw;
            cw.bx = bx.data(); cw.by = by.data(); cw.bz = bz.data();
            cw.dx = f.x.data(); cw.dy = f.y.data(); cw.dz = f.z.data();
            cw.ux = u.x.data(); cw.uy = u.y.data(); cw.uz = u.z.data();
            cw.ras2vox = levelR2V; cw.nx = lnx; cw.ny = lny; cw.nz = lnz;
            cw.ox = tmp.x.data(); cw.oy = tmp.y.data(); cw.oz = tmp.z.data();
            TinyParallel::parallelFor(0, N, cw, 1024);
          }
          std::swap(u.x, tmp.x); std::swap(u.y, tmp.y); std::swap(u.z, tmp.z);
          if (totalSigma > 0.0) {
            smoothFieldInto(u, totalSigma, gscratch.data());
            zeroFieldBoundary(u);
          }
        }

        trace.push_back(metricVal);
        ++itRun;
        lastCost = metricVal;
        double conv = R_PosInf;
        if (monitorOn) {
          conv = monitor.add(energy0);
          if (conv < convergenceThreshold) {
            stopCode = 1;
          } else if (maxnSide[0] <= 1e-12 && maxnSide[1] <= 1e-12 && totalSigma <= 0.0) {
            stopCode = 2;
          }
        }
        lastConv = conv;
        if (verbose) {
          MaxNormReducer a(curF.x.data(), curF.y.data(), curF.z.data());
          TinyParallel::parallelReduce(0, N, a, 4096);
          MaxNormReducer b(curM.x.data(), curM.y.data(), curM.z.data());
          TinyParallel::parallelReduce(0, N, b, 4096);
          Rcpp::Rcout << "  it " << std::setw(5) << (it + 1) << ": cost = " << std::fixed
                      << std::setprecision(6) << metricVal << "  field_max (fixed/moving) = "
                      << std::setprecision(2) << std::sqrt(a.maxn2) << "/" << std::sqrt(b.maxn2) << " mm";
          if (monitorOn && std::isfinite(conv)) {
            Rcpp::Rcout << "  convergence = " << std::scientific << std::setprecision(3) << conv << std::defaultfloat;
          }
          Rcpp::Rcout << std::endl;
        }
        Rcpp::checkUserInterrupt();
        if (stopCode != 0) break;
      }
      iterationsRun[lev] = itRun;
      if (verbose && monitorOn) {
        Rcpp::Rcout << "  => ";
        if (stopCode == 1) {
          Rcpp::Rcout << "converged (convergence " << std::scientific << std::setprecision(3) << lastConv
                      << " < threshold " << convergenceThreshold << std::defaultfloat << ") after " << itRun
                      << " of " << maxIter << " iterations";
        } else if (stopCode == 2) {
          Rcpp::Rcout << "force vanished after " << itRun << " of " << maxIter << " iterations";
        } else {
          Rcpp::Rcout << "max iterations reached (" << maxIter << ")";
        }
        Rcpp::Rcout << "  cost=" << std::fixed << std::setprecision(6) << lastCost << std::endl;
      }
    } // maxIter > 0

    dispF = curF; dispM = curM; haveField = true;
  }

  // ---- full grid: phi_F^-1 (inverted), composed with uM ----
  const std::size_t Nf = (std::size_t)fnx * fny * fnz;
  const std::size_t szf = (std::size_t)fnx * (std::size_t)fny;
  Field full; full.nx = fnx; full.ny = fny; full.nz = fnz; full.vox2ras = fV2R;
  full.alloc(fnx, fny, fnz);
  ravereg_apply::InvertOptions iopt;
  iopt.maxIterations = 20;
  if (haveField) {
    Field fullF; fullF.nx = fnx; fullF.ny = fny; fullF.nz = fnz; fullF.vox2ras = fV2R;
    fullF.alloc(fnx, fny, fnz);
    resampleField(dispF, fullF);
    dispF = Field();
    std::vector<double> inv(3 * Nf, 0.0);
    ravereg_apply::invertDisplacementField<reg_real>(
      fullF.x.data(), fullF.y.data(), fullF.z.data(), fnx, fny, fnz, fV2R,
      inv.data(), inv.data() + Nf, inv.data() + 2 * Nf, iopt);
    fullF = Field();
    std::vector<reg_real> ix(Nf), iy(Nf), iz(Nf);
    for (std::size_t i = 0; i < Nf; ++i) {
      ix[i] = (reg_real)inv[i]; iy[i] = (reg_real)inv[Nf + i]; iz[i] = (reg_real)inv[2 * Nf + i];
    }
    std::vector<double>().swap(inv);
    Field fullM; fullM.nx = fnx; fullM.ny = fny; fullM.nz = fnz; fullM.vox2ras = fV2R;
    fullM.alloc(fnx, fny, fnz);
    resampleField(dispM, fullM);
    dispM = Field();
    std::vector<double> bx(Nf), by(Nf), bz(Nf);
    {
      BaseGridWorker bw;
      bw.v2r = fV2R; bw.nx = fnx; bw.sz = szf;
      bw.bx = bx.data(); bw.by = by.data(); bw.bz = bz.data();
      TinyParallel::parallelFor(0, Nf, bw, 1024);
    }
    ComposeWorker<reg_real> cw;
    cw.bx = bx.data(); cw.by = by.data(); cw.bz = bz.data();
    cw.dx = ix.data(); cw.dy = iy.data(); cw.dz = iz.data();
    cw.ux = fullM.x.data(); cw.uy = fullM.y.data(); cw.uz = fullM.z.data();
    cw.ras2vox = fV2R.inverse(); cw.nx = fnx; cw.ny = fny; cw.nz = fnz;
    cw.ox = full.x.data(); cw.oy = full.y.data(); cw.oz = full.z.data();
    TinyParallel::parallelFor(0, Nf, cw, 1024);
  }

  // ---- outputs, as the one-sided path ----
  std::vector<std::vector<reg_real> > movStore(K);
  std::vector<RegImage3> movFulls; movFulls.reserve(K);
  for (int kk = 0; kk < K; ++kk) {
    movStore[kk].assign(movingVecs[kk].begin(), movingVecs[kk].end());
    movFulls.push_back(RegImage3(movStore[kk].data(), mnx, mny, mnz, mV2R));
  }
  std::vector<Rcpp::NumericVector> warpedCh(K);
  std::vector<double*> outPtrs(K);
  std::vector<int> codes(K);
  for (int kk = 0; kk < K; ++kk) {
    warpedCh[kk] = Rcpp::NumericVector(Nf);
    outPtrs[kk] = &warpedCh[kk][0];
    codes[kk] = (kk < interpCodes.size()) ? interpCodes[kk] : 1;
  }
  {
    OutputWarpWorker ow;
    ow.fV2R = fV2R; ow.aff = aff;
    ow.dispx = full.x.data(); ow.dispy = full.y.data(); ow.dispz = full.z.data();
    ow.fnx = fnx; ow.szf = szf; ow.movFulls = movFulls.data(); ow.codes = codes.data();
    ow.K = K; ow.outs = outPtrs.data();
    TinyParallel::parallelFor(0, Nf, ow, 256);
  }
  Rcpp::List images(K);
  for (int kk = 0; kk < K; ++kk) {
    warpedCh[kk].attr("dim") = Rcpp::IntegerVector::create(fnx, fny, fnz);
    images[kk] = warpedCh[kk];
  }
  Rcpp::NumericVector fwd(Nf * 3);
  std::copy(full.x.begin(), full.x.end(), fwd.begin());
  std::copy(full.y.begin(), full.y.end(), fwd.begin() + Nf);
  std::copy(full.z.begin(), full.z.end(), fwd.begin() + 2 * Nf);
  fwd.attr("dim") = Rcpp::IntegerVector::create(fnx, fny, fnz, 3);
  Rcpp::NumericVector inverseField(static_cast<R_xlen_t>(3 * Nf));
  if (haveField) {
    double* v = &inverseField[0];
    const ravereg_apply::InvertResult ires = ravereg_apply::invertDisplacementField<reg_real>(
      full.x.data(), full.y.data(), full.z.data(), fnx, fny, fnz, fV2R, v, v + Nf, v + 2 * Nf, iopt);
    if (verbose) {
      Rcpp::Rcout << "[SyN, symmetric] inverse field: " << ires.iterations << " fixed-point iterations"
                  << ", residual mean = " << std::scientific << std::setprecision(2) << ires.meanResidual
                  << ", max = " << ires.maxResidual << " voxel" << std::defaultfloat << std::endl;
    }
  }
  inverseField.attr("dim") = Rcpp::IntegerVector::create(fnx, fny, fnz, 3);
  return Rcpp::List::create(
    Rcpp::Named("forward_field") = fwd,
    Rcpp::Named("inverse_field") = inverseField,
    Rcpp::Named("image") = warpedCh[0],
    Rcpp::Named("images") = images,
    Rcpp::Named("metric_trace") = Rcpp::wrap(trace),
    Rcpp::Named("iterations") = iterationsRun);
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List register_syn_cpp(const Rcpp::List& fixedList,
                            const Rcpp::IntegerVector& fixedDim,
                            const Rcpp::NumericMatrix& fixedVox2Ras,
                            const Rcpp::List& movingList,
                            const Rcpp::IntegerVector& movingDim,
                            const Rcpp::NumericMatrix& movingVox2Ras,
                            const Rcpp::NumericMatrix& affineTransform,
                            const Rcpp::CharacterVector& metricNames,
                            const Rcpp::NumericVector& weights,
                            const Rcpp::IntegerVector& shrinkFactors,
                            const Rcpp::NumericVector& smoothingSigmas,
                            const Rcpp::IntegerVector& iterations,
                            const double gradStep,
                            const Rcpp::NumericVector& flowSigma,
                            const double totalSigma,
                            const int ccRadius,
                            const Rcpp::IntegerVector& interpCodes,
                            const Rcpp::Nullable<Rcpp::NumericVector>& fixedMask,
                            const Rcpp::Nullable<Rcpp::NumericVector>& movingMask,
                            const Rcpp::NumericMatrix& fixedPoints,
                            const Rcpp::NumericMatrix& movingPoints,
                            const double pointsWeight,
                            const double convergenceThreshold,
                            const int convergenceWindow,
                            const int miBins = 32,
                            const double miSampling = 0.25,
                            const bool normalizeChannels = false,
                            const int seed = 1,
                            const bool symmetric = false,
                            const bool verbose = false) {
  // Optional per-level convergence stop (ConvergenceMonitor): a window of at
  // least 2 energies and a finite threshold; window 0 (or a threshold that is
  // not finite) runs every iteration
  if (convergenceWindow != 0 && convergenceWindow < 2) {
    Rcpp::stop("register_syn_cpp: `convergenceWindow` must be 0 (no stop) or at least 2");
  }
  const bool monitorOn = convergenceWindow >= 2 && R_finite(convergenceThreshold);

  const int fnx = fixedDim[0], fny = fixedDim[1], fnz = fixedDim[2];
  const int mnx = movingDim[0], mny = movingDim[1], mnz = movingDim[2];
  const Matrix4d fV2R = asMatrix4(fixedVox2Ras);
  const Matrix4d mV2R = asMatrix4(movingVox2Ras);
  const Matrix4d aff  = asMatrix4(affineTransform);
  const Matrix3d affR = aff.topLeftCorner<3, 3>();
  const Matrix3d affRT = affR.transpose();

  const int K = fixedList.size();

  // Keep channel data alive for the whole call; record per-channel metric/weight.
  std::vector<Rcpp::NumericVector> fixedVecs, movingVecs;
  fixedVecs.reserve(K); movingVecs.reserve(K);
  // metric of each channel: 0 mean squares, 1 local CC, 2 Mattes MI
  std::vector<int> mcode(K);
  std::vector<double> w(K);
  bool anyMI = false;
  for (int k = 0; k < K; ++k) {
    fixedVecs.push_back(Rcpp::as<Rcpp::NumericVector>(fixedList[k]));
    movingVecs.push_back(Rcpp::as<Rcpp::NumericVector>(movingList[k]));
    const std::string mn = Rcpp::as<std::string>(metricNames[k]);
    mcode[k] = (mn == "meansquares") ? 0 : ((mn == "mattes") ? 2 : 1);
    if (mcode[k] == 2) anyMI = true;
    w[k] = weights[k];
  }
  if (anyMI && miBins < 2 * kMIPadding + 4) {
    Rcpp::stop("register_syn_cpp: `miBins` must be at least 8");
  }

  // Optional masks. The moving mask (full moving grid) is sampled at warped
  // points; the fixed mask (full fixed grid) is shrunk + dilated per level.
  std::vector<reg_real> movMaskStore;
  RegImage3 movMaskImg;
  const RegImage3* movMaskPtr = nullptr;
  if (movingMask.isNotNull()) {
    Rcpp::NumericVector mm(movingMask.get());
    movMaskStore.assign(mm.begin(), mm.end());
    movMaskImg = RegImage3(movMaskStore.data(), mnx, mny, mnz, mV2R);
    movMaskPtr = &movMaskImg;
  }
  std::vector<reg_real> fixedMaskFull;
  const reg_real* fixedMaskData = nullptr;
  if (fixedMask.isNotNull()) {
    Rcpp::NumericVector fm(fixedMask.get());
    fixedMaskFull.assign(fm.begin(), fm.end());
    fixedMaskData = fixedMaskFull.data();
  }

  // Optional cortical landmark correspondences (N x 3 RAS each; N=0 => none).
  // fixedPoints are target/fixed-RAS anchors; movingPoints the corresponding
  // source/moving-RAS targets (row i = the same cortical vertex). They add a
  // weighted force pulling warp(fixedPoints[i]) toward movingPoints[i].
  const int nP = fixedPoints.nrow();
  std::vector<Vector3d> fixPts, movPts;
  if (nP > 0 && pointsWeight > 0.0) {
    fixPts.reserve(nP); movPts.reserve(nP);
    for (int p = 0; p < nP; ++p) {
      fixPts.emplace_back(fixedPoints(p, 0), fixedPoints(p, 1), fixedPoints(p, 2));
      movPts.emplace_back(movingPoints(p, 0), movingPoints(p, 1), movingPoints(p, 2));
    }
  }

  // symmetric SyN: its own loop (no landmarks)
  if (symmetric) {
    if (!fixPts.empty()) {
      Rcpp::stop("register_syn_cpp: landmarks are not available with the symmetric deformable stage");
    }
    return registerSynSymmetric(
      fixedVecs, movingVecs, fnx, fny, fnz, mnx, mny, mnz, fV2R, mV2R, aff, mcode, w,
      shrinkFactors, smoothingSigmas, iterations, gradStep, flowSigma, totalSigma, ccRadius,
      interpCodes, fixedMaskData, movMaskPtr, movMaskStore, miBins, miSampling,
      normalizeChannels, seed, monitorOn, convergenceThreshold, convergenceWindow, verbose);
  }

  const int nLevels = shrinkFactors.size();
  std::vector<double> trace;
  Rcpp::IntegerVector iterationsRun(nLevels);   // per level; 0 when skipped

  Field disp;        // forward field on the current fixed level
  bool haveField = false;

  for (int lev = 0; lev < nLevels; ++lev) {
    const int shrink = shrinkFactors[lev];
    const double sigma = (lev < smoothingSigmas.size()) ? smoothingSigmas[lev] : 0.0;
    const int maxIter = (lev < iterations.size()) ? iterations[lev] : 0;
    // per-level fluid (flow) regularization sigma; clamp to the last entry
    const double flowSigmaLev = (flowSigma.size() == 0) ? 0.0
        : (lev < flowSigma.size() ? flowSigma[lev] : flowSigma[flowSigma.size() - 1]);

    // Per-channel smoothed images + level structures (all channels share the
    // grid). The smoothed arrays are borrowed by FixedLevel/RegImage3 below, so
    // they must outlive them -> keep in level-scoped storage.
    std::vector<std::vector<reg_real> > sFixedStore(K), sMovingStore(K);
    std::vector<FixedLevel> flvls; flvls.reserve(K);
    std::vector<RegImage3> movs; movs.reserve(K);
    for (int k = 0; k < K; ++k) {
      sFixedStore[k]  = gaussianSmooth(&fixedVecs[k][0], fnx, fny, fnz, sigma);
      sMovingStore[k] = gaussianSmooth(&movingVecs[k][0], mnx, mny, mnz, sigma);
      flvls.push_back(buildFixedLevel(sFixedStore[k].data(), fnx, fny, fnz, fV2R, shrink));
      movs.push_back(RegImage3(sMovingStore[k].data(), mnx, mny, mnz, mV2R));
    }

    const FixedLevel& flvl0 = flvls[0];
    const int lnx = flvl0.nx, lny = flvl0.ny, lnz = flvl0.nz;
    const std::size_t N = flvl0.size();
    const std::size_t szl = (std::size_t)lnx * (std::size_t)lny;
    const Matrix4d levelV2R = flvl0.vox2ras;

    // (re)build the shared field at this level
    Field cur; cur.nx = lnx; cur.ny = lny; cur.nz = lnz; cur.vox2ras = levelV2R;
    cur.alloc(lnx, lny, lnz);
    if (haveField) {
      resampleField(disp, cur);        // cur is already allocated/zeroed
    }

    const double levelSpacing = std::cbrt(std::abs(
        levelV2R.topLeftCorner<3, 3>().determinant()));

    if (verbose) {
      Rcpp::Rcout << "[SyN] level " << (lev+1) << "/" << nLevels
                  << " (shrink=" << shrink << ", sigma=" << std::fixed
                  << std::setprecision(1) << sigma
                  << ", flow=" << std::setprecision(1) << flowSigmaLev
                  << ", channels=" << K << ")"
                  << ": max " << maxIter << " iterations" << std::endl;
    }

    // Everything below is only used by the iteration loop; skip the work (and
    // the large allocations) entirely on levels with no deformable iterations
    // (e.g. the default finest level, syn_iterations = c(40, 20, 0)).
    if (maxIter > 0) {
      const double winN = std::pow(2.0 * ccRadius + 1.0, 3);

      // Per-channel local-CC precomputed fixed sums (F is constant per level).
      std::vector<std::vector<double> > sumF(K), sumFF(K);
      std::vector<double> ccTau(K, 0.0);
      for (int k = 0; k < K; ++k) {
        if (mcode[k] != 1) continue;
        const std::vector<reg_real>& Fk = flvls[k].data;
        sumF[k] = boxSum(Fk, lnx, lny, lnz, ccRadius);
        std::vector<double> FF(N);
        for (std::size_t i = 0; i < N; ++i) FF[i] = (double)Fk[i] * (double)Fk[i];
        sumFF[k] = boxSum(FF, lnx, lny, lnz, ccRadius);
        // global fixed variance -> threshold below which a window is "flat"
        double gm = 0.0;
        for (std::size_t i = 0; i < N; ++i) gm += Fk[i];
        gm /= (double) N;
        double gv = 0.0;
        for (std::size_t i = 0; i < N; ++i) { const double d = Fk[i] - gm; gv += d * d; }
        gv /= (double) N;
        ccTau[k] = 1e-3 * gv * winN;
      }

      // Precompute the per-level base RAS grid (constant across iterations).
      std::vector<double> bx(N), by(N), bz(N);
      {
        BaseGridWorker bw;
        bw.v2r = levelV2R; bw.nx = lnx; bw.sz = szl;
        bw.bx = bx.data(); bw.by = by.data(); bw.bz = bz.data();
        TinyParallel::parallelFor(0, N, bw, 1024);
      }

      // Per-level fixed mask (shrunk to 0/1) and its ccRadius dilation. maskLev
      // gates which voxels receive force; warpMask = (boxSum(maskLev) > 0) gates
      // the costly warp (a voxel is warped iff some force voxel's window needs
      // its value), so background voxels are skipped while every force voxel's
      // CC window still sees correct warped values.
      std::vector<reg_real> maskLev, warpMask;
      if (fixedMaskData) {
        FixedLevel ml = buildFixedLevel(fixedMaskData, fnx, fny, fnz, fV2R, shrink);
        maskLev.resize(N);
        for (std::size_t i = 0; i < N; ++i)
          maskLev[i] = (ml.data[i] > reg_real(0.5)) ? reg_real(1) : reg_real(0);
        std::vector<double> md = boxSum(maskLev, lnx, lny, lnz, ccRadius);
        warpMask.resize(N);
        for (std::size_t i = 0; i < N; ++i)
          warpMask[i] = (md[i] > 0.0) ? reg_real(1) : reg_real(0);
      }
      const reg_real* maskLevPtr  = maskLev.empty()  ? nullptr : maskLev.data();
      const reg_real* warpMaskPtr = warpMask.empty() ? nullptr : warpMask.data();

      // Hoisted per-iteration scratch (reused every iteration / channel).
      std::vector<double> Mw(N), gFx(N), gFy(N), gFz(N);
      std::vector<double> sumM(N), sumMM(N), sumFM(N);   // local-CC window sums
      std::vector<double> scrM(N), scrMM(N), scrFM(N);   // their box-sum intermediates
      std::vector<reg_real> gscratch(N);                 // Gaussian ping-pong scratch
      Field force; force.nx = lnx; force.ny = lny; force.nz = lnz;
      force.vox2ras = levelV2R; force.alloc(lnx, lny, lnz);
      // block partial sums of the cost (fixed blocks, combined in order)
      const std::size_t nCostBlocks = (N + kCostBlock - 1) / kCostBlock;
      std::vector<double> blockAcc(nCostBlocks), blockSse(nCostBlocks);
      std::vector<long> blockCnt(nCostBlocks);

      // Normalized channel combination: each channel's force goes to its own
      // field first (allocated only when used)
      Field chForce;
      std::vector<double> normBlocks;
      if (normalizeChannels) {
        chForce.nx = lnx; chForce.ny = lny; chForce.nz = lnz;
        chForce.vox2ras = levelV2R; chForce.alloc(lnx, lny, lnz);
        normBlocks.assign(nCostBlocks, 0.0);
      }

      // Mattes MI channels: one seeded random sample of the level voxels (inside
      // the fixed mask), shared by the channels, and the intensity range of
      // each channel (the fixed level image inside the fixed mask, the smoothed
      // moving image inside the moving mask)
      std::vector<std::size_t> miSamples;
      std::vector<uint8_t> miValid;
      std::vector<MIRange> miRange(K);
      std::vector<double> miJoint, miLog;
      if (anyMI) {
        miValid.assign(N, 0);
        std::mt19937 rng((unsigned int)seed + 7919u * (unsigned int)(lev + 1));
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        for (std::size_t i = 0; i < N; ++i) {
          if (maskLevPtr && maskLevPtr[i] == reg_real(0)) continue;
          if (miSampling >= 1.0 || unif(rng) < miSampling) miSamples.push_back(i);
        }
        for (int k = 0; k < K; ++k) {
          if (mcode[k] != 2) continue;
          const std::vector<reg_real>& Fk = flvls[k].data;
          double f0 = R_PosInf, f1 = R_NegInf;
          for (std::size_t i = 0; i < N; ++i) {
            if (maskLevPtr && maskLevPtr[i] == reg_real(0)) continue;
            const double v = Fk[i];
            if (v < f0) f0 = v;
            if (v > f1) f1 = v;
          }
          const std::vector<reg_real>& Mk = sMovingStore[k];
          double m0 = R_PosInf, m1 = R_NegInf;
          for (std::size_t i = 0; i < Mk.size(); ++i) {
            if (movMaskPtr && movMaskStore[i] <= reg_real(0.5)) continue;
            const double v = Mk[i];
            if (v < m0) m0 = v;
            if (v > m1) m1 = v;
          }
          if (!std::isfinite(f0)) { f0 = 0.0; f1 = 1.0; }
          if (!std::isfinite(m0)) { m0 = 0.0; m1 = 1.0; }
          miRange[k].set(f0, f1, m0, m1, miBins);
        }
      }

      // Noise floor of the normalized combination: a channel force below 1e-6
      // of the channel's own scale (its moving-gradient energy times a metric
      // factor: the fixed spread for mean squares, 2 / (spread * window) for
      // local CC, 1 / Δm for MI) is rounding noise, such as a channel that is
      // already aligned, which the normalization would blow up to its weight
      std::vector<double> chScale(K, 0.0);
      if (normalizeChannels) {
        for (int k = 0; k < K; ++k) {
          if (mcode[k] == 2) { chScale[k] = 1.0 / miRange[k].mbin; continue; }
          const std::vector<reg_real>& Fk = flvls[k].data;
          double s1 = 0.0, s2 = 0.0; std::size_t n = 0;
          for (std::size_t i = 0; i < N; ++i) {
            if (maskLevPtr && maskLevPtr[i] == reg_real(0)) continue;
            const double v = Fk[i]; s1 += v; s2 += v * v; ++n;
          }
          double sd = 0.0;
          if (n > 1) {
            const double m = s1 / (double)n;
            sd = std::sqrt(std::max(0.0, s2 / (double)n - m * m));
          }
          chScale[k] = (mcode[k] == 0) ? sd : ((sd > 0.0) ? 2.0 / (sd * winN) : 0.0);
        }
      }

      // a fresh convergence monitor at every level, as in ITK
      ConvergenceMonitor monitor(monitorOn ? convergenceWindow : 0);
      int itRun = 0, stopCode = 0;      // stopCode: 1 converged, 2 force vanished
      double lastCost = 0.0, lastConv = R_PosInf;

      for (int it = 0; it < maxIter; ++it) {
        force.zero();
        double metricVal = 0.0;
        double energy0 = 0.0;           // unweighted cost of the first channel

        double totalMag = 0.0;          // normalized combination: sum of |d_k|
        for (int k = 0; k < K; ++k) {
          const std::vector<reg_real>& Fk = flvls[k].data;
          RegImage3& movk = movs[k];
          const double wk = w[k];
          // where the channel force goes, and with which weight: the shared
          // force field with wk, or (normalized) the channel field with 1
          reg_real* tx = force.x.data();
          reg_real* ty = force.y.data();
          reg_real* tz = force.z.data();
          double wForce = wk;
          if (normalizeChannels) {
            chForce.zero();
            tx = chForce.x.data(); ty = chForce.y.data(); tz = chForce.z.data();
            wForce = 1.0;
          }

          // ---- warp moving channel k + spatial gradient (fixed-space) ----
          {
            WarpWorker ww;
            ww.mov = &movk; ww.aff = aff; ww.affRT = affRT;
            ww.movMask = movMaskPtr; ww.warpMask = warpMaskPtr;
            ww.bx = bx.data(); ww.by = by.data(); ww.bz = bz.data();
            ww.dispx = cur.x.data(); ww.dispy = cur.y.data(); ww.dispz = cur.z.data();
            ww.Mw = Mw.data(); ww.gFx = gFx.data(); ww.gFy = gFy.data(); ww.gFz = gFz.data();
            ww.valid = (mcode[k] == 2) ? miValid.data() : nullptr;
            TinyParallel::parallelFor(0, N, ww, 256);
          }

          // ---- channel similarity force, accumulated weighted ----
          if (mcode[k] == 1) {
            fusedLocalCCSums(Mw.data(), Fk.data(), lnx, lny, lnz, ccRadius,
                             sumM, sumMM, sumFM, scrM, scrMM, scrFM);

            CCForceWorker cc;
            cc.F = Fk.data(); cc.Mw = Mw.data();
            cc.sumF = sumF[k].data(); cc.sumFF = sumFF[k].data();
            cc.sumM = sumM.data(); cc.sumMM = sumMM.data(); cc.sumFM = sumFM.data();
            cc.gFx = gFx.data(); cc.gFy = gFy.data(); cc.gFz = gFz.data();
            cc.fx = tx; cc.fy = ty; cc.fz = tz;
            cc.winN = winN; cc.tau = ccTau[k]; cc.wk = wForce;
            cc.mask = maskLevPtr; cc.N = N;
            cc.blockAcc = blockAcc.data(); cc.blockCnt = blockCnt.data();
            TinyParallel::parallelFor(0, nCostBlocks, cc, 1);
            double accCC = 0.0; long cnt = 0;
            for (std::size_t b = 0; b < nCostBlocks; ++b) {
              accCC += blockAcc[b]; cnt += blockCnt[b];
            }
            const double ek = (cnt > 0) ? -(accCC / cnt) : 0.0;   // cost = -localCC
            if (k == 0) energy0 = ek;
            metricVal += wk * ek;
          } else if (mcode[k] == 2) {
            const double mi = mattesForce(Fk.data(), Mw.data(), miValid.data(),
                                          gFx.data(), gFy.data(), gFz.data(),
                                          miSamples, miRange[k], wForce, tx, ty, tz,
                                          miJoint, miLog);
            const double ek = -mi;                                   // cost = -MI
            if (k == 0) energy0 = ek;
            metricVal += wk * ek;
          } else {
            MSForceWorker ms;
            ms.F = Fk.data(); ms.Mw = Mw.data();
            ms.gFx = gFx.data(); ms.gFy = gFy.data(); ms.gFz = gFz.data();
            ms.fx = tx; ms.fy = ty; ms.fz = tz;
            ms.wk = wForce; ms.mask = maskLevPtr; ms.N = N;
            ms.blockSse = blockSse.data();
            TinyParallel::parallelFor(0, nCostBlocks, ms, 1);
            double sse = 0.0;
            for (std::size_t b = 0; b < nCostBlocks; ++b) sse += blockSse[b];
            const double ek = sse / N;
            if (k == 0) energy0 = ek;
            metricVal += wk * ek;
          }

          // ---- normalized combination: add w_k d_k / |d_k| ----
          if (normalizeChannels) {
            const double nk = std::sqrt(blockSumSquares(
              chForce.x.data(), chForce.y.data(), chForce.z.data(), N, normBlocks));
            const double gk = std::sqrt(blockSumSquares(
              gFx.data(), gFy.data(), gFz.data(), N, normBlocks));
            if (nk > 1e-30 && nk > 1e-6 * gk * chScale[k]) {
              totalMag += nk;
              AxpyFieldWorker aw;
              aw.dx = force.x.data(); aw.dy = force.y.data(); aw.dz = force.z.data();
              aw.sx = chForce.x.data(); aw.sy = chForce.y.data(); aw.sz = chForce.z.data();
              aw.a = wk / nk;
              TinyParallel::parallelFor(0, N, aw, 4096);
            }
          }
        } // channels
        // as ITK, rescale by the summed magnitudes (keeps the image force on the
        // scale of the landmark term below; the step cap sets the update size)
        if (normalizeChannels && totalMag > 0.0) {
          ScaleFieldWorker sw;
          sw.fx = force.x.data(); sw.fy = force.y.data(); sw.fz = force.z.data();
          sw.a = totalMag;
          TinyParallel::parallelFor(0, N, sw, 4096);
        }

        // ---- cortical landmark force (surface-guided), splatted onto force ----
        // For each correspondence, pull warp(fixed anchor) toward its moving
        // target: gather the current displacement at the anchor, form the
        // residual error, convert to displacement-space force (affR^T e), and
        // scatter it onto the field. The fluid smoothing below then turns these
        // sparse pulls into a smooth deformation, exactly like the image force.
        if (!fixPts.empty()) {
          const Matrix4d s2v = levelV2R.inverse();
          double surfCost = 0.0;
          for (int p = 0; p < nP; ++p) {
            const Vector3d& xf = fixPts[p];
            Vector4d vh = s2v * Vector4d(xf[0], xf[1], xf[2], 1.0);
            const double cx = vh[0], cy = vh[1], cz = vh[2];
            double dx, dy, dz;
            const bool ok =
              trilinearSample(cur.x.data(), lnx, lny, lnz, cx, cy, cz, dx) &&
              trilinearSample(cur.y.data(), lnx, lny, lnz, cx, cy, cz, dy) &&
              trilinearSample(cur.z.data(), lnx, lny, lnz, cx, cy, cz, dz);
            if (!ok) continue;                       // anchor outside this grid
            Vector4d xdh(xf[0] + dx, xf[1] + dy, xf[2] + dz, 1.0);
            Vector3d q = (aff * xdh).head<3>();
            Vector3d e = movPts[p] - q;
            surfCost += 0.5 * e.squaredNorm();
            Vector3d fvec = affRT * e;               // -d(1/2|e|^2)/d disp
            splatTrilinear(force, cx, cy, cz, pointsWeight, fvec);
          }
          metricVal += pointsWeight * (surfCost / (double)nP);
        }

        // ---- regularize update (fluid) ----
        smoothFieldInto(force, flowSigmaLev, gscratch.data());

        // cap the step to grad_step * voxel spacing
        MaxNormReducer mn(force.x.data(), force.y.data(), force.z.data());
        TinyParallel::parallelReduce(0, N, mn, 4096);
        const double maxn = std::sqrt(mn.maxn2);
        const double scale = (maxn > 1e-12) ? (gradStep * levelSpacing / maxn) : 0.0;

        {
          FieldUpdateWorker uw;
          uw.cx = cur.x.data(); uw.cy = cur.y.data(); uw.cz = cur.z.data();
          uw.fx = force.x.data(); uw.fy = force.y.data(); uw.fz = force.z.data();
          uw.scale = scale;
          TinyParallel::parallelFor(0, N, uw, 4096);
        }

        // ---- regularize total field (elastic) ----
        smoothFieldInto(cur, totalSigma, gscratch.data());

        trace.push_back(metricVal);
        ++itRun;
        lastCost = metricVal;

        // ---- convergence, after the update as in ITK ----
        // (the cost was measured before this iteration's update)
        double conv = R_PosInf;
        if (monitorOn) {
          conv = monitor.add(energy0);
          if (conv < convergenceThreshold) {
            stopCode = 1;
          } else if (maxn <= 1e-12 && totalSigma <= 0.0) {
            stopCode = 2;     // no update was applied: the field cannot change anymore
          }
        }
        lastConv = conv;

        if (verbose) {
          MaxNormReducer fm(cur.x.data(), cur.y.data(), cur.z.data());
          TinyParallel::parallelReduce(0, N, fm, 4096);
          Rcpp::Rcout << "  it " << std::setw(5) << (it+1)
                      << ": cost = " << std::fixed << std::setprecision(6) << metricVal
                      << "  field_max = " << std::setprecision(2)
                      << std::sqrt(fm.maxn2) << " mm";
          if (monitorOn && std::isfinite(conv)) {
            Rcpp::Rcout << "  convergence = " << std::scientific << std::setprecision(3)
                        << conv << std::defaultfloat;
          }
          Rcpp::Rcout << std::endl;
        }

        Rcpp::checkUserInterrupt();
        if (stopCode != 0) break;
      }
      iterationsRun[lev] = itRun;

      if (verbose && monitorOn) {
        Rcpp::Rcout << "  => ";
        if (stopCode == 1) {
          Rcpp::Rcout << "converged (convergence " << std::scientific << std::setprecision(3)
                      << lastConv << " < threshold " << convergenceThreshold
                      << std::defaultfloat << ") after " << itRun << " of " << maxIter
                      << " iterations";
        } else if (stopCode == 2) {
          Rcpp::Rcout << "force vanished after " << itRun << " of " << maxIter << " iterations";
        } else {
          Rcpp::Rcout << "max iterations reached (" << maxIter << ")";
        }
        Rcpp::Rcout << "  cost=" << std::fixed << std::setprecision(6) << lastCost
                    << std::endl;
      }
    } // maxIter > 0

    disp = cur; haveField = true;
  }

  // ---- warp every moving channel onto the full-resolution fixed grid ----
  Field full; full.nx = fnx; full.ny = fny; full.nz = fnz; full.vox2ras = fV2R;
  full.alloc(fnx, fny, fnz);
  if (haveField) resampleField(disp, full);

  // RegImage3 stores reg_real; convert each (unsmoothed) moving channel.
  std::vector<std::vector<reg_real> > movStore(K);
  std::vector<RegImage3> movFulls; movFulls.reserve(K);
  for (int kk = 0; kk < K; ++kk) {
    movStore[kk].assign(movingVecs[kk].begin(), movingVecs[kk].end());
    movFulls.push_back(RegImage3(movStore[kk].data(), mnx, mny, mnz, mV2R));
  }

  const std::size_t Nf = full.size();
  const std::size_t szf = (std::size_t)fnx * (std::size_t)fny;

  std::vector<Rcpp::NumericVector> warpedCh(K);
  std::vector<double*> outPtrs(K);
  std::vector<int> codes(K);
  for (int kk = 0; kk < K; ++kk) {
    warpedCh[kk] = Rcpp::NumericVector(Nf);
    outPtrs[kk] = &warpedCh[kk][0];
    codes[kk] = (kk < interpCodes.size()) ? interpCodes[kk] : 1;
  }

  {
    OutputWarpWorker ow;
    ow.fV2R = fV2R; ow.aff = aff;
    ow.dispx = full.x.data(); ow.dispy = full.y.data(); ow.dispz = full.z.data();
    ow.fnx = fnx; ow.szf = szf;
    ow.movFulls = movFulls.data();
    ow.codes = codes.data();
    ow.K = K;
    ow.outs = outPtrs.data();
    TinyParallel::parallelFor(0, Nf, ow, 256);
  }

  Rcpp::List images(K);
  for (int kk = 0; kk < K; ++kk) {
    warpedCh[kk].attr("dim") = Rcpp::IntegerVector::create(fnx, fny, fnz);
    images[kk] = warpedCh[kk];
  }

  // pack the forward / inverse fields as (nx, ny, nz, 3) arrays
  auto packField = [&](Field& f) {
    Rcpp::NumericVector v(f.size() * 3);
    std::copy(f.x.begin(), f.x.end(), v.begin());
    std::copy(f.y.begin(), f.y.end(), v.begin() + f.size());
    std::copy(f.z.begin(), f.z.end(), v.begin() + 2 * f.size());
    v.attr("dim") = Rcpp::IntegerVector::create(f.nx, f.ny, f.nz, 3);
    return v;
  };

  // ---- inverse field: fixed-point inversion of the final forward field ----
  // (nx, ny, nz, 3) in double, on the same full fixed grid as `full`.
  Rcpp::NumericVector inverseField(static_cast<R_xlen_t>(3 * Nf));
  if (haveField) {
    double* v = &inverseField[0];
    ravereg_apply::InvertOptions iopt;
    // at most 20 sweeps, as the inversion inside ANTs' SyN: nodes at the border
    // of an unmasked field often have no inverse and never meet the tolerances,
    // while the interior is far below them after 20 sweeps
    iopt.maxIterations = 20;
    const ravereg_apply::InvertResult ires = ravereg_apply::invertDisplacementField<reg_real>(
      full.x.data(), full.y.data(), full.z.data(), fnx, fny, fnz, fV2R,
      v, v + Nf, v + 2 * Nf, iopt);
    if (verbose) {
      Rcpp::Rcout << "[SyN] inverse field: " << ires.iterations << " fixed-point iterations"
                  << ", residual mean = " << std::scientific << std::setprecision(2)
                  << ires.meanResidual << ", max = " << ires.maxResidual
                  << " voxel" << std::defaultfloat << std::endl;
    }
  }
  inverseField.attr("dim") = Rcpp::IntegerVector::create(fnx, fny, fnz, 3);

  return Rcpp::List::create(
    Rcpp::Named("forward_field") = packField(full),
    Rcpp::Named("inverse_field") = inverseField,
    Rcpp::Named("image") = warpedCh[0],     // primary channel (back-compat)
    Rcpp::Named("images") = images,         // all channels, per-channel interpolation
    Rcpp::Named("metric_trace") = Rcpp::wrap(trace),
    Rcpp::Named("iterations") = iterationsRun);   // iterations run per level
}

// The convergence values of the deformable-stage stop for the energy trace of
// one level (one value per iteration, +Inf before the window is full), through
// the same ConvergenceMonitor as register_syn_cpp. Internal, for the tests.
// [[Rcpp::export]]
Rcpp::NumericVector syn_convergence_values_cpp(const Rcpp::NumericVector& energy,
                                               const int window) {
  if (window < 2) {
    Rcpp::stop("syn_convergence_values_cpp: `window` must be at least 2");
  }
  ConvergenceMonitor monitor(window);
  Rcpp::NumericVector out(energy.size());
  for (R_xlen_t k = 0; k < energy.size(); ++k) out[k] = monitor.add(energy[k]);
  return out;
}
