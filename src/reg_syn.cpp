// Deformable registration: greedy symmetric diffeomorphic ("SyN"-style).
//
// After a linear (affine) stage, a dense RAS displacement field on the FIXED
// grid captures the residual nonlinear deformation. Each iteration:
//   1. warp the moving image to the fixed grid through affine o (id + disp);
//   2. form a similarity force (local normalized CC, or mean squares);
//   3. Gaussian-smooth the update (flow_sigma), the fluid-like regularizer;
//   4. take a capped step (grad_step), updating the forward field and, with the
//      opposite update, the inverse field (symmetric);
//   5. optionally Gaussian-smooth the total field (total_sigma), elastic.
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
// field is bit-identical regardless of thread count; only the scalar cost trace
// (a parallel sum) may differ in its last ulps.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "reg_core.h"

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

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      // Outside the dilated fixed mask there is nothing to drive and no force
      // voxel needs this value in its window: skip the (dominant) sampling work.
      if (warpMask && warpMask[o] == reg_real(0)) {
        Mw[o] = 0.0; gFx[o] = 0.0; gFy[o] = 0.0; gFz[o] = 0.0;
        continue;
      }
      Vector4d xdh(bx[o] + dispx[o], by[o] + dispy[o], bz[o] + dispz[o], 1.0);
      Vector3d q = (aff * xdh).head<3>();
      double mval; Vector3d gM;
      mov->sampleRas(q, mval, &gM);
      // A point mapping outside the moving mask is treated as background (no
      // value, no gradient) so it contributes neither force nor structure.
      if (movMask) {
        double mv; movMask->sampleRas(q, mv);
        if (mv <= 0.5) { Mw[o] = 0.0; gFx[o] = 0.0; gFy[o] = 0.0; gFz[o] = 0.0; continue; }
      }
      Mw[o] = mval;
      Vector3d gF = affRT * gM;
      gFx[o] = gF[0]; gFy[o] = gF[1]; gFz[o] = gF[2];
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

// Local normalized-CC force, accumulated (weighted) into the shared force
// field; the scalar local-CC sum is reduced for the cost trace. Each index i is
// written by one thread (disjoint ranges), so the force is bit-identical to the
// serial result; only accCC/cnt are order-dependent.
struct CCForceReducer : public TinyParallel::Worker {
  const reg_real* F;
  const double *Mw, *sumF, *sumFF, *sumM, *sumMM, *sumFM;
  const double *gFx, *gFy, *gFz;
  reg_real *fx, *fy, *fz;
  double winN, tau, wk;
  const reg_real* mask = nullptr;   // optional fixed mask: skip force where 0
  double accCC; long cnt;

  CCForceReducer(const reg_real* F_, const double* Mw_,
                 const double* sumF_, const double* sumFF_, const double* sumM_,
                 const double* sumMM_, const double* sumFM_,
                 const double* gFx_, const double* gFy_, const double* gFz_,
                 reg_real* fx_, reg_real* fy_, reg_real* fz_,
                 double winN_, double tau_, double wk_)
    : F(F_), Mw(Mw_), sumF(sumF_), sumFF(sumFF_), sumM(sumM_), sumMM(sumMM_),
      sumFM(sumFM_), gFx(gFx_), gFy(gFy_), gFz(gFz_),
      fx(fx_), fy(fy_), fz(fz_), winN(winN_), tau(tau_), wk(wk_),
      accCC(0.0), cnt(0) {}

  CCForceReducer(CCForceReducer& o, TinyParallel::Split)
    : F(o.F), Mw(o.Mw), sumF(o.sumF), sumFF(o.sumFF), sumM(o.sumM),
      sumMM(o.sumMM), sumFM(o.sumFM), gFx(o.gFx), gFy(o.gFy), gFz(o.gFz),
      fx(o.fx), fy(o.fy), fz(o.fz), winN(o.winN), tau(o.tau), wk(o.wk),
      mask(o.mask), accCC(0.0), cnt(0) {}

  void operator()(std::size_t begin, std::size_t end) override {
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
  }

  void join(const CCForceReducer& rhs) { accCC += rhs.accCC; cnt += rhs.cnt; }
};

// Mean-squares force, accumulated (weighted) into the shared force field; the
// SSE is reduced for the cost trace.
struct MSForceReducer : public TinyParallel::Worker {
  const reg_real* F;
  const double *Mw, *gFx, *gFy, *gFz;
  reg_real *fx, *fy, *fz;
  double wk;
  const reg_real* mask = nullptr;   // optional fixed mask: skip force where 0
  double sse;

  MSForceReducer(const reg_real* F_, const double* Mw_,
                 const double* gFx_, const double* gFy_, const double* gFz_,
                 reg_real* fx_, reg_real* fy_, reg_real* fz_, double wk_)
    : F(F_), Mw(Mw_), gFx(gFx_), gFy(gFy_), gFz(gFz_),
      fx(fx_), fy(fy_), fz(fz_), wk(wk_), sse(0.0) {}

  MSForceReducer(MSForceReducer& o, TinyParallel::Split)
    : F(o.F), Mw(o.Mw), gFx(o.gFx), gFy(o.gFy), gFz(o.gFz),
      fx(o.fx), fy(o.fy), fz(o.fz), wk(o.wk), mask(o.mask), sse(0.0) {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      if (mask && mask[i] == reg_real(0)) continue;   // outside fixed mask
      const double d = (double)F[i] - Mw[i];
      fx[i] += wk * d * gFx[i];
      fy[i] += wk * d * gFy[i];
      fz[i] += wk * d * gFz[i];
      sse += d * d;
    }
  }

  void join(const MSForceReducer& rhs) { sse += rhs.sse; }
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

// cur += scale * force ; curInv -= scale * force  (symmetric update).
struct FieldUpdateWorker : public TinyParallel::Worker {
  reg_real *cx, *cy, *cz, *ix, *iy, *iz;
  const reg_real *fx, *fy, *fz;
  double scale;
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      const double sx = scale * fx[i], sy = scale * fy[i], sz = scale * fz[i];
      cx[i] += sx; cy[i] += sy; cz[i] += sz;
      ix[i] -= sx; iy[i] -= sy; iz[i] -= sz;
    }
  }
};

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
                            const bool verbose = false) {
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
  std::vector<char> useCC(K);
  std::vector<double> w(K);
  for (int k = 0; k < K; ++k) {
    fixedVecs.push_back(Rcpp::as<Rcpp::NumericVector>(fixedList[k]));
    movingVecs.push_back(Rcpp::as<Rcpp::NumericVector>(movingList[k]));
    const std::string mn = Rcpp::as<std::string>(metricNames[k]);
    useCC[k] = (mn != "meansquares") ? 1 : 0;
    w[k] = weights[k];
  }

  // Optional masks. The moving mask (full moving grid) is sampled at warped
  // points; the fixed mask (full fixed grid) is shrunk + dilated per level.
  std::vector<reg_real> movMaskStore;
  RegImage3 movMaskImg;
  const RegImage3* movMaskPtr = nullptr;
  if (movingMask.isNotNull()) {
    Rcpp::NumericVector mm(movingMask);
    movMaskStore.assign(mm.begin(), mm.end());
    movMaskImg = RegImage3(movMaskStore.data(), mnx, mny, mnz, mV2R);
    movMaskPtr = &movMaskImg;
  }
  std::vector<reg_real> fixedMaskFull;
  const reg_real* fixedMaskData = nullptr;
  if (fixedMask.isNotNull()) {
    Rcpp::NumericVector fm(fixedMask);
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

  const int nLevels = shrinkFactors.size();
  std::vector<double> trace;

  Field disp;        // forward field on the current fixed level
  Field dispInv;     // approximate inverse field
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
    Field curInv = cur;
    if (haveField) {
      resampleField(disp, cur);        // cur is already allocated/zeroed
      resampleField(dispInv, curInv);
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
        if (!useCC[k]) continue;
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

      for (int it = 0; it < maxIter; ++it) {
        force.zero();
        double metricVal = 0.0;

        for (int k = 0; k < K; ++k) {
          const std::vector<reg_real>& Fk = flvls[k].data;
          RegImage3& movk = movs[k];
          const double wk = w[k];

          // ---- warp moving channel k + spatial gradient (fixed-space) ----
          {
            WarpWorker ww;
            ww.mov = &movk; ww.aff = aff; ww.affRT = affRT;
            ww.movMask = movMaskPtr; ww.warpMask = warpMaskPtr;
            ww.bx = bx.data(); ww.by = by.data(); ww.bz = bz.data();
            ww.dispx = cur.x.data(); ww.dispy = cur.y.data(); ww.dispz = cur.z.data();
            ww.Mw = Mw.data(); ww.gFx = gFx.data(); ww.gFy = gFy.data(); ww.gFz = gFz.data();
            TinyParallel::parallelFor(0, N, ww, 256);
          }

          // ---- channel similarity force, accumulated weighted ----
          if (useCC[k]) {
            fusedLocalCCSums(Mw.data(), Fk.data(), lnx, lny, lnz, ccRadius,
                             sumM, sumMM, sumFM, scrM, scrMM, scrFM);

            CCForceReducer cc(
              Fk.data(), Mw.data(), sumF[k].data(), sumFF[k].data(),
              sumM.data(), sumMM.data(), sumFM.data(),
              gFx.data(), gFy.data(), gFz.data(),
              force.x.data(), force.y.data(), force.z.data(),
              winN, ccTau[k], wk);
            cc.mask = maskLevPtr;
            TinyParallel::parallelReduce(0, N, cc, 1024);
            metricVal += wk * ((cc.cnt > 0) ? -(cc.accCC / cc.cnt) : 0.0);   // cost = -localCC
          } else {
            MSForceReducer ms(
              Fk.data(), Mw.data(), gFx.data(), gFy.data(), gFz.data(),
              force.x.data(), force.y.data(), force.z.data(), wk);
            ms.mask = maskLevPtr;
            TinyParallel::parallelReduce(0, N, ms, 1024);
            metricVal += wk * (ms.sse / N);
          }
        } // channels

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
          uw.ix = curInv.x.data(); uw.iy = curInv.y.data(); uw.iz = curInv.z.data();
          uw.fx = force.x.data(); uw.fy = force.y.data(); uw.fz = force.z.data();
          uw.scale = scale;
          TinyParallel::parallelFor(0, N, uw, 4096);
        }

        // ---- regularize total field (elastic) ----
        smoothFieldInto(cur, totalSigma, gscratch.data());
        smoothFieldInto(curInv, totalSigma, gscratch.data());

        trace.push_back(metricVal);

        if (verbose) {
          MaxNormReducer fm(cur.x.data(), cur.y.data(), cur.z.data());
          TinyParallel::parallelReduce(0, N, fm, 4096);
          Rcpp::Rcout << "  it " << std::setw(5) << (it+1)
                      << ": cost = " << std::fixed << std::setprecision(6) << metricVal
                      << "  field_max = " << std::setprecision(2)
                      << std::sqrt(fm.maxn2) << " mm" << std::endl;
        }

        Rcpp::checkUserInterrupt();
      }
    } // maxIter > 0

    disp = cur; dispInv = curInv; haveField = true;
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

  Field fullInv; fullInv.nx = fnx; fullInv.ny = fny; fullInv.nz = fnz;
  fullInv.vox2ras = fV2R; fullInv.alloc(fnx, fny, fnz);
  if (haveField) resampleField(dispInv, fullInv);

  return Rcpp::List::create(
    Rcpp::Named("forward_field") = packField(full),
    Rcpp::Named("inverse_field") = packField(fullInv),
    Rcpp::Named("image") = warpedCh[0],     // primary channel (back-compat)
    Rcpp::Named("images") = images,         // all channels, per-channel interpolation
    Rcpp::Named("metric_trace") = Rcpp::wrap(trace));
}
