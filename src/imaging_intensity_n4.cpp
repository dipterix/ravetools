// Native N4 bias-field correction (Tustison et al. 2010, IEEE TMI 29(6)) and
// the two preprocessing helpers that the 'abpN4' recipe of ANTsR / ANTsPy
// (0.3.x) wraps around it: histogram-quantile intensity truncation and the
// morphological mask clean-up of `get_mask()`.
//
// This is a clean reimplementation from the published algorithms (N3 histogram
// sharpening: Sled, Zijdenbos & Evans 1998; multilevel B-spline approximation:
// Lee, Wolberg & Shin 1997). ITK / ANTs were consulted only as behavioral
// references for defaults and conventions (padding rule, shrink sampling,
// kernel centering, convergence measure); no code was copied.
//
// Conventions
//   * volumes are column-major (x fastest), 0-indexed; all arithmetic double
//   * the B-spline lattice lives on a *padded* grid whose extent is a whole
//     number of spline spans (the ANTs padding rule); the padding is never
//     materialized, only its geometry (lower bound, padded size) is tracked
//   * the fit runs on the sub-sampled ("shrunk") grid, the final bias field is
//     evaluated on the full-resolution grid, both through the same parametric
//     domain, so fit and reconstruction are geometrically consistent
//   * every TinyParallel section writes disjoint output ranges, and partial
//     sums are reduced in a fixed order, so results are bit-identical for any
//     thread count
//   * no R API is touched inside worker threads; Rcpp::checkUserInterrupt()
//     runs once per (serial) N4 iteration

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include "TinyParallel.h"

namespace raven4 {

typedef std::size_t usize;

// ---------------------------------------------------------------------------
// 1. Intensity truncation quantiles (ANTs ImageMath TruncateImageIntensity)
//
// Histogram of the finite voxels > 0 with `bins` equal-width bins spanning
// their [min, max]; the quantile is linearly interpolated inside the bin that
// crosses the requested probability, walking up from the bottom for p < 0.5
// and down from the top otherwise (ITK Histogram::Quantile semantics).
// ---------------------------------------------------------------------------
static double histogram_quantile(const std::vector<double>& H, double total,
                                 double minv, double maxv, double p)
{
  const int bins = static_cast<int>(H.size());
  const double interval = (maxv - minv) / static_cast<double>(bins);
  auto bin_lo = [&](int j) { return minv + static_cast<double>(j) * interval; };
  auto bin_hi = [&](int j) {
    return (j >= bins - 1) ? maxv : minv + static_cast<double>(j + 1) * interval;
  };
  if (p < 0.5) {
    int n = 0;
    double cum = 0.0, p_n = 0.0, p_prev = 0.0, f = 0.0;
    do {
      f = H[n];
      cum += f;
      p_prev = p_n;
      p_n = cum / total;
      ++n;
    } while (n < bins && p_n < p);
    const int j = n - 1;
    const double prop = f / total;
    const double lo = bin_lo(j), hi = bin_hi(j);
    if (!(prop > 0.0)) return lo;
    return lo + ((p - p_prev) / prop) * (hi - lo);
  } else {
    int n = bins - 1, m = 0;
    double cum = 0.0, p_n = 1.0, p_prev = 1.0, f = 0.0;
    do {
      f = H[n];
      cum += f;
      p_prev = p_n;
      p_n = 1.0 - cum / total;
      --n;
      ++m;
    } while (m < bins && p_n > p);
    const int j = n + 1;
    const double prop = f / total;
    const double lo = bin_lo(j), hi = bin_hi(j);
    if (!(prop > 0.0)) return hi;
    return hi - ((p_prev - p) / prop) * (hi - lo);
  }
}

// ---------------------------------------------------------------------------
// 2. Binary morphology on a 0/1 volume (ITK conventions)
//
// The structuring element is the ITK "ball": all integer offsets d with
// |d|^2 <= (radius + 0.5)^2. Erosion treats the outside of the volume as
// foreground (border voxels are not eroded by the image edge), dilation treats
// it as background. Both parallelize over z slices (independent outputs).
// ---------------------------------------------------------------------------
struct BallOffsets {
  std::vector<int> dx, dy, dz;
};

static BallOffsets ball_offsets(int radius)
{
  BallOffsets b;
  const double r2 = (static_cast<double>(radius) + 0.5) * (static_cast<double>(radius) + 0.5);
  for (int z = -radius; z <= radius; ++z)
    for (int y = -radius; y <= radius; ++y)
      for (int x = -radius; x <= radius; ++x) {
        if (x == 0 && y == 0 && z == 0) continue;
        const double d2 = static_cast<double>(x * x + y * y + z * z);
        if (d2 <= r2) { b.dx.push_back(x); b.dy.push_back(y); b.dz.push_back(z); }
      }
  return b;
}

struct MorphWorker : public TinyParallel::Worker {
  const unsigned char* in;
  unsigned char* out;
  int nx, ny, nz;
  const BallOffsets* off;
  bool erode;

  MorphWorker(const unsigned char* in_, unsigned char* out_, int nx_, int ny_, int nz_,
              const BallOffsets* off_, bool erode_)
    : in(in_), out(out_), nx(nx_), ny(ny_), nz(nz_), off(off_), erode(erode_) {}

  void operator()(usize begin, usize end) override {
    const usize sy = static_cast<usize>(nx);
    const usize sz = static_cast<usize>(nx) * static_cast<usize>(ny);
    const usize no = off->dx.size();
    for (usize zz = begin; zz < end; ++zz) {
      const int z = static_cast<int>(zz);
      for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
          const usize v = static_cast<usize>(x) + sy * static_cast<usize>(y) + sz * zz;
          if (erode) {
            if (!in[v]) { out[v] = 0; continue; }
            unsigned char keep = 1;
            for (usize o = 0; o < no; ++o) {
              const int xx = x + off->dx[o], yy = y + off->dy[o], zq = z + off->dz[o];
              if (xx < 0 || yy < 0 || zq < 0 || xx >= nx || yy >= ny || zq >= nz) continue;
              if (!in[static_cast<usize>(xx) + sy * static_cast<usize>(yy) + sz * static_cast<usize>(zq)]) {
                keep = 0; break;
              }
            }
            out[v] = keep;
          } else {
            if (in[v]) { out[v] = 1; continue; }
            unsigned char hit = 0;
            for (usize o = 0; o < no; ++o) {
              const int xx = x + off->dx[o], yy = y + off->dy[o], zq = z + off->dz[o];
              if (xx < 0 || yy < 0 || zq < 0 || xx >= nx || yy >= ny || zq >= nz) continue;
              if (in[static_cast<usize>(xx) + sy * static_cast<usize>(yy) + sz * static_cast<usize>(zq)]) {
                hit = 1; break;
              }
            }
            out[v] = hit;
          }
        }
      }
    }
  }
};

static std::vector<unsigned char> morph(const std::vector<unsigned char>& in,
                                        int nx, int ny, int nz, int radius, bool erode)
{
  std::vector<unsigned char> out(in.size(), 0);
  if (radius <= 0) { out = in; return out; }
  const BallOffsets off = ball_offsets(radius);
  MorphWorker w(in.data(), out.data(), nx, ny, nz, &off, erode);
  TinyParallel::parallelFor(0, static_cast<usize>(nz), w, 1);
  return out;
}

// Face-connected (6-neighborhood) components of the voxels equal to `target`.
// labels: 0 for voxels not in the set, 1..k otherwise; sizes[k-1] = voxel count.
static int label_components(const std::vector<unsigned char>& m, int nx, int ny, int nz,
                            unsigned char target, std::vector<int>& labels,
                            std::vector<usize>& sizes)
{
  const usize n = m.size();
  const usize sy = static_cast<usize>(nx);
  const usize sz = static_cast<usize>(nx) * static_cast<usize>(ny);
  labels.assign(n, 0);
  sizes.clear();
  std::vector<usize> stack;
  int k = 0;
  for (usize seed = 0; seed < n; ++seed) {
    if (m[seed] != target || labels[seed] != 0) continue;
    ++k;
    usize count = 0;
    labels[seed] = k;
    stack.clear();
    stack.push_back(seed);
    while (!stack.empty()) {
      const usize v = stack.back();
      stack.pop_back();
      ++count;
      const usize z = v / sz;
      const usize r = v - z * sz;
      const usize y = r / sy;
      const usize x = r - y * sy;
      usize nb[6];
      int nn = 0;
      if (x > 0) nb[nn++] = v - 1;
      if (x + 1 < static_cast<usize>(nx)) nb[nn++] = v + 1;
      if (y > 0) nb[nn++] = v - sy;
      if (y + 1 < static_cast<usize>(ny)) nb[nn++] = v + sy;
      if (z > 0) nb[nn++] = v - sz;
      if (z + 1 < static_cast<usize>(nz)) nb[nn++] = v + sz;
      for (int i = 0; i < nn; ++i) {
        const usize u = nb[i];
        if (m[u] == target && labels[u] == 0) { labels[u] = k; stack.push_back(u); }
      }
    }
    sizes.push_back(count);
  }
  return k;
}

// Keep the largest foreground component(s); components with fewer than
// `min_size` voxels are discarded first (ANTs GetLargestComponent).
static std::vector<unsigned char> largest_component(const std::vector<unsigned char>& m,
                                                    int nx, int ny, int nz, usize min_size)
{
  std::vector<int> labels;
  std::vector<usize> sizes;
  label_components(m, nx, ny, nz, 1, labels, sizes);
  usize best = 0;
  for (usize s : sizes) if (s >= min_size && s > best) best = s;
  std::vector<unsigned char> out(m.size(), 0);
  if (best == 0) return out;
  for (usize v = 0; v < m.size(); ++v) {
    const int l = labels[v];
    if (l > 0 && sizes[static_cast<usize>(l - 1)] == best) out[v] = 1;
  }
  return out;
}

// Fill holes: every face-connected background component except the largest
// one becomes foreground (ANTs FillHoles with hole parameter 2).
static std::vector<unsigned char> fill_holes(const std::vector<unsigned char>& m,
                                             int nx, int ny, int nz)
{
  std::vector<int> labels;
  std::vector<usize> sizes;
  const int k = label_components(m, nx, ny, nz, 0, labels, sizes);
  std::vector<unsigned char> out(m);
  if (k <= 1) return out;
  int largest = 1;
  for (int i = 2; i <= k; ++i) {
    if (sizes[static_cast<usize>(i - 1)] > sizes[static_cast<usize>(largest - 1)]) largest = i;
  }
  for (usize v = 0; v < m.size(); ++v) {
    const int l = labels[v];
    if (l > 0 && l != largest) out[v] = 1;
  }
  return out;
}

static std::vector<unsigned char> logical_to_bytes(const Rcpp::LogicalVector& mask)
{
  std::vector<unsigned char> m(static_cast<usize>(mask.size()), 0);
  for (R_xlen_t i = 0; i < mask.size(); ++i) m[static_cast<usize>(i)] = (mask[i] == TRUE) ? 1 : 0;
  return m;
}

static Rcpp::LogicalVector bytes_to_logical(const std::vector<unsigned char>& m)
{
  Rcpp::LogicalVector out(static_cast<R_xlen_t>(m.size()));
  for (usize i = 0; i < m.size(); ++i) out[static_cast<R_xlen_t>(i)] = m[i] ? TRUE : FALSE;
  return out;
}

static void check_dims(const Rcpp::IntegerVector& dims, R_xlen_t n, const char* what)
{
  if (dims.size() != 3) Rcpp::stop("%s: `dims` must have length 3.", what);
  for (int d = 0; d < 3; ++d) {
    if (dims[d] == NA_INTEGER || dims[d] < 1) Rcpp::stop("%s: invalid dimensions.", what);
  }
  const double prod = static_cast<double>(dims[0]) * static_cast<double>(dims[1]) *
                      static_cast<double>(dims[2]);
  if (prod != static_cast<double>(n)) Rcpp::stop("%s: `dims` do not match the data length.", what);
}

// ---------------------------------------------------------------------------
// 3. Uniform B-spline machinery
//
// A lattice of control points phi[i, j, k] (column-major) spans a parametric
// domain [0, S_d] per axis with S_d = ncp_d - order spans. Control point i is
// centered at u = i - (order - 1) / 2 and the basis of the (order + 1)
// control points floor(u) + a, a = 0..order, is K_order(u - floor(u) - a +
// (order - 1) / 2) with K the centered cardinal B-spline of that degree.
// ---------------------------------------------------------------------------
static inline double bspline_kernel(int order, double v)
{
  const double a = std::fabs(v);
  switch (order) {
  case 1:
    return (a < 1.0) ? 1.0 - a : 0.0;
  case 2:
    if (a < 0.5) return 0.75 - a * a;
    if (a < 1.5) return (9.0 - 12.0 * a + 4.0 * a * a) * 0.125;
    return 0.0;
  default:  // 3
    if (a < 1.0) return (4.0 - 6.0 * a * a + 3.0 * a * a * a) / 6.0;
    if (a < 2.0) return (8.0 - 12.0 * a + 6.0 * a * a - a * a * a) / 6.0;
    return 0.0;
  }
}

// Per-axis sampling table: for each sample along the axis, the first control
// index of its support, the (order + 1) basis weights and their sum of squares.
struct AxisTable {
  int n = 0;
  int order = 3;
  std::vector<int> start;
  std::vector<double> w;    // n * (order + 1)
  std::vector<double> sq;   // n
  const double* weights(int i) const { return &w[static_cast<usize>(i) * static_cast<usize>(order + 1)]; }
};

// `u` holds the parametric coordinate of every sample in [0, spans].
static AxisTable make_axis_table(int order, int spans, const std::vector<double>& u)
{
  AxisTable t;
  t.n = static_cast<int>(u.size());
  t.order = order;
  t.start.assign(u.size(), 0);
  t.w.assign(u.size() * static_cast<usize>(order + 1), 0.0);
  t.sq.assign(u.size(), 0.0);
  const double half = 0.5 * static_cast<double>(order - 1);
  for (usize i = 0; i < u.size(); ++i) {
    int span = static_cast<int>(std::floor(u[i]));
    if (span > spans - 1) span = spans - 1;
    if (span < 0) span = 0;
    const double frac = u[i] - static_cast<double>(span);
    t.start[i] = span;
    double s2 = 0.0;
    for (int a = 0; a <= order; ++a) {
      const double b = bspline_kernel(order, frac - static_cast<double>(a) + half);
      t.w[i * static_cast<usize>(order + 1) + static_cast<usize>(a)] = b;
      s2 += b * b;
    }
    t.sq[i] = s2;
  }
  return t;
}

struct Lattice {
  int nx = 0, ny = 0, nz = 0;
  std::vector<double> v;
  usize size() const { return static_cast<usize>(nx) * static_cast<usize>(ny) * static_cast<usize>(nz); }
  void alloc(int a, int b, int c) { nx = a; ny = b; nz = c; v.assign(size(), 0.0); }
};

// Evaluate the lattice on the grid described by three axis tables, one z slice
// per work item. The lattice is collapsed one dimension at a time (z, then y,
// then x), which costs about (order + 1) multiply-adds per output voxel.
struct ReconWorker : public TinyParallel::Worker {
  const Lattice* lat;
  const AxisTable *tx, *ty, *tz;
  double* out;

  ReconWorker(const Lattice* lat_, const AxisTable* tx_, const AxisTable* ty_,
              const AxisTable* tz_, double* out_)
    : lat(lat_), tx(tx_), ty(ty_), tz(tz_), out(out_) {}

  void operator()(usize begin, usize end) override {
    const int ncx = lat->nx, ncy = lat->ny;
    const int order = tx->order;
    const usize nx = static_cast<usize>(tx->n), ny = static_cast<usize>(ty->n);
    std::vector<double> L2(static_cast<usize>(ncx) * static_cast<usize>(ncy));
    std::vector<double> L1(static_cast<usize>(ncx));
    const double* phi = lat->v.data();
    for (usize z = begin; z < end; ++z) {
      const int sz = tz->start[z];
      const double* wz = tz->weights(static_cast<int>(z));
      for (int jy = 0; jy < ncy; ++jy) {
        for (int ix = 0; ix < ncx; ++ix) {
          double acc = 0.0;
          for (int c = 0; c <= order; ++c) {
            acc += wz[c] * phi[static_cast<usize>(ix) + static_cast<usize>(ncx) *
                   (static_cast<usize>(jy) + static_cast<usize>(ncy) * static_cast<usize>(sz + c))];
          }
          L2[static_cast<usize>(ix) + static_cast<usize>(ncx) * static_cast<usize>(jy)] = acc;
        }
      }
      for (usize y = 0; y < ny; ++y) {
        const int sy = ty->start[y];
        const double* wy = ty->weights(static_cast<int>(y));
        for (int ix = 0; ix < ncx; ++ix) {
          double acc = 0.0;
          for (int b = 0; b <= order; ++b) {
            acc += wy[b] * L2[static_cast<usize>(ix) + static_cast<usize>(ncx) * static_cast<usize>(sy + b)];
          }
          L1[static_cast<usize>(ix)] = acc;
        }
        double* row = out + nx * (y + ny * z);
        for (usize x = 0; x < nx; ++x) {
          const int sx = tx->start[x];
          const double* wx = tx->weights(static_cast<int>(x));
          double acc = 0.0;
          for (int a = 0; a <= order; ++a) acc += wx[a] * L1[static_cast<usize>(sx + a)];
          row[x] = acc;
        }
      }
    }
  }
};

static void reconstruct(const Lattice& lat, const AxisTable& tx, const AxisTable& ty,
                        const AxisTable& tz, double* out)
{
  ReconWorker w(&lat, &tx, &ty, &tz, out);
  TinyParallel::parallelFor(0, static_cast<usize>(tz.n), w, 1);
}

// Lattice refinement (Lane-Riesenfeld subdivision): doubling the number of
// spans, the refined control points 2m + off are weighted sums of the coarse
// points m + j with the binomial masks of the uniform B-spline of that degree.
// The refined lattice reproduces exactly the same function.
static Lattice refine_lattice(const Lattice& in, int order)
{
  double c[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
  switch (order) {
  case 1:  c[0][0] = 1.0;   c[1][0] = 0.5;   c[1][1] = 0.5; break;
  case 2:  c[0][0] = 0.75;  c[0][1] = 0.25;  c[1][0] = 0.25;  c[1][1] = 0.75; break;
  default: c[0][0] = 0.5;   c[0][1] = 0.5;   c[1][0] = 0.125; c[1][1] = 0.75; c[1][2] = 0.125; break;
  }
  Lattice out;
  out.alloc(2 * in.nx - order, 2 * in.ny - order, 2 * in.nz - order);
  for (int kz = 0; kz < out.nz; ++kz) {
    const int mz = kz / 2, oz = kz % 2;
    for (int ky = 0; ky < out.ny; ++ky) {
      const int my = ky / 2, oy = ky % 2;
      for (int kx = 0; kx < out.nx; ++kx) {
        const int mx = kx / 2, ox = kx % 2;
        double sum = 0.0;
        for (int jz = 0; jz <= order; ++jz) {
          if (mz + jz >= in.nz) break;
          const double cz = c[oz][jz];
          if (cz == 0.0) continue;
          for (int jy = 0; jy <= order; ++jy) {
            if (my + jy >= in.ny) break;
            const double cyz = c[oy][jy] * cz;
            if (cyz == 0.0) continue;
            for (int jx = 0; jx <= order; ++jx) {
              if (mx + jx >= in.nx) break;
              const double cx = c[ox][jx];
              if (cx == 0.0) continue;
              sum += cx * cyz * in.v[static_cast<usize>(mx + jx) + static_cast<usize>(in.nx) *
                     (static_cast<usize>(my + jy) + static_cast<usize>(in.ny) * static_cast<usize>(mz + jz))];
            }
          }
        }
        out.v[static_cast<usize>(kx) + static_cast<usize>(out.nx) *
              (static_cast<usize>(ky) + static_cast<usize>(out.ny) * static_cast<usize>(kz))] = sum;
      }
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// 4. Scattered-data B-spline approximation (Lee, Wolberg & Shin 1997)
//
// For every included voxel p with value z_p and confidence w_p, each control
// point c in its support receives
//     omega_c += w_p B_c(p)^2
//     delta_c += w_p B_c(p)^3 z_p / sum_{c'} B_{c'}(p)^2
// and phi_c = delta_c / omega_c (0 where nothing contributes). The volume is
// split into a fixed number of z-chunks, each accumulating its own partial
// lattices in a fixed voxel order; the partials are then summed in chunk
// order, so the result does not depend on the thread count.
// ---------------------------------------------------------------------------
struct FitChunkWorker : public TinyParallel::Worker {
  const unsigned char* inc;
  const double* res;
  const double* wt;
  int mx, my;
  const AxisTable *tx, *ty, *tz;
  int ncx, ncy;
  usize Ls;
  const std::vector<int>* chunk_z0;
  const std::vector<int>* chunk_z1;
  double* omega_part;
  double* delta_part;

  FitChunkWorker(const unsigned char* inc_, const double* res_, const double* wt_,
                 int mx_, int my_, const AxisTable* tx_, const AxisTable* ty_,
                 const AxisTable* tz_, int ncx_, int ncy_, usize Ls_,
                 const std::vector<int>* z0_, const std::vector<int>* z1_,
                 double* om_, double* de_)
    : inc(inc_), res(res_), wt(wt_), mx(mx_), my(my_), tx(tx_), ty(ty_), tz(tz_),
      ncx(ncx_), ncy(ncy_), Ls(Ls_), chunk_z0(z0_), chunk_z1(z1_),
      omega_part(om_), delta_part(de_) {}

  void operator()(usize begin, usize end) override {
    const int order = tx->order;
    const usize smx = static_cast<usize>(mx), smy = static_cast<usize>(my);
    for (usize g = begin; g < end; ++g) {
      double* om = omega_part + g * Ls;
      double* de = delta_part + g * Ls;
      const int z0 = (*chunk_z0)[g], z1 = (*chunk_z1)[g];
      for (int z = z0; z < z1; ++z) {
        const int sz = tz->start[z];
        const double* wz = tz->weights(z);
        const double sqz = tz->sq[z];
        for (int y = 0; y < my; ++y) {
          const int sy = ty->start[y];
          const double* wy = ty->weights(y);
          const double sqyz = ty->sq[y] * sqz;
          const usize rowbase = smx * (static_cast<usize>(y) + smy * static_cast<usize>(z));
          for (int x = 0; x < mx; ++x) {
            const usize v = rowbase + static_cast<usize>(x);
            if (!inc[v]) continue;
            const int sx = tx->start[x];
            const double* wx = tx->weights(x);
            const double w2sum = tx->sq[x] * sqyz;
            const double wgt = wt[v];
            const double cc = res[v] * wgt / w2sum;
            for (int c = 0; c <= order; ++c) {
              const double Bz = wz[c];
              const usize kb = static_cast<usize>(sz + c);
              for (int b = 0; b <= order; ++b) {
                const double Byz = wy[b] * Bz;
                const usize row = static_cast<usize>(ncx) * (static_cast<usize>(sy + b) + static_cast<usize>(ncy) * kb) +
                                  static_cast<usize>(sx);
                for (int a = 0; a <= order; ++a) {
                  const double B = wx[a] * Byz;
                  const double B2 = B * B;
                  om[row + static_cast<usize>(a)] += wgt * B2;
                  de[row + static_cast<usize>(a)] += cc * B2 * B;
                }
              }
            }
          }
        }
      }
    }
  }
};

// Reduce the chunk partials (in chunk order) and add the resulting lattice
// increment to `phi`.
struct FitReduceWorker : public TinyParallel::Worker {
  const double* omega_part;
  const double* delta_part;
  usize G, Ls;
  double* phi;

  FitReduceWorker(const double* om_, const double* de_, usize G_, usize Ls_, double* phi_)
    : omega_part(om_), delta_part(de_), G(G_), Ls(Ls_), phi(phi_) {}

  void operator()(usize begin, usize end) override {
    for (usize e = begin; e < end; ++e) {
      double om = 0.0, de = 0.0;
      for (usize g = 0; g < G; ++g) {
        om += omega_part[g * Ls + e];
        de += delta_part[g * Ls + e];
      }
      double p = 0.0;
      if (om > 0.0) {
        p = de / om;
        if (!std::isfinite(p)) p = 0.0;
      }
      phi[e] += p;
    }
  }
};

// ---------------------------------------------------------------------------
// 5. Histogram sharpening (Sled et al. 1998, as parameterized in N4)
// ---------------------------------------------------------------------------
typedef std::complex<double> cplx;

// In-place iterative radix-2 FFT; `a.size()` must be a power of two. The
// inverse is normalized by 1/N.
static void fft_inplace(std::vector<cplx>& a, bool inverse)
{
  const usize n = a.size();
  for (usize i = 1, j = 0; i < n; ++i) {
    usize bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  const double pi = std::acos(-1.0);
  for (usize len = 2; len <= n; len <<= 1) {
    const usize half = len >> 1;
    const double ang = (inverse ? 2.0 : -2.0) * pi / static_cast<double>(len);
    for (usize i = 0; i < n; i += len) {
      for (usize j = 0; j < half; ++j) {
        const cplx w(std::cos(ang * static_cast<double>(j)), std::sin(ang * static_cast<double>(j)));
        const cplx u = a[i + j];
        const cplx v = a[i + j + half] * w;
        a[i + j] = u + v;
        a[i + j + half] = u - v;
      }
    }
  }
  if (inverse) {
    const double s = 1.0 / static_cast<double>(n);
    for (usize i = 0; i < n; ++i) a[i] *= s;
  }
}

// Map the log intensities of the included voxels through the sharpened
// histogram expectation E(u | v): deconvolve the (triangular-binned) histogram
// of `vals` by a Gaussian of the given FWHM with a Wiener filter, then compute
// the conditional expectation of the true intensity given the observed one.
// Writes the sharpened value of every included voxel into `out` (same index
// list). If the intensities are degenerate the values pass through unchanged.
static void sharpen_histogram(const double* vals, const std::vector<usize>& idx,
                              int bins, double fwhm, double noise, double* out)
{
  double minv = std::numeric_limits<double>::infinity();
  double maxv = -std::numeric_limits<double>::infinity();
  for (usize v : idx) {
    const double p = vals[v];
    if (p < minv) minv = p;
    if (p > maxv) maxv = p;
  }
  const double slope = (maxv - minv) / static_cast<double>(bins - 1);
  if (!(slope > 0.0) || !std::isfinite(slope)) {
    for (usize v : idx) out[v] = vals[v];
    return;
  }

  std::vector<double> H(static_cast<usize>(bins), 0.0);
  for (usize v : idx) {
    const double cidx = (vals[v] - minv) / slope;
    int b = static_cast<int>(std::floor(cidx));
    if (b < 0) b = 0;
    if (b > bins - 1) b = bins - 1;
    const double offset = cidx - static_cast<double>(b);
    if (offset == 0.0) {
      H[static_cast<usize>(b)] += 1.0;
    } else if (b < bins - 1) {
      H[static_cast<usize>(b)] += 1.0 - offset;
      H[static_cast<usize>(b + 1)] += offset;
    } else {
      H[static_cast<usize>(b)] += 1.0;
    }
  }

  // zero-pad to twice the smallest power of two that holds the histogram
  usize padded = 1;
  while (padded < static_cast<usize>(bins)) padded <<= 1;
  padded <<= 1;
  const usize hoff = (padded - static_cast<usize>(bins)) / 2;

  std::vector<cplx> V(padded, cplx(0.0, 0.0));
  for (int b = 0; b < bins; ++b) V[hoff + static_cast<usize>(b)] = cplx(H[static_cast<usize>(b)], 0.0);
  std::vector<cplx> Vf(V);
  fft_inplace(Vf, false);

  // Gaussian blur kernel (periodic, centered at 0) and its spectrum
  const double ln2 = std::log(2.0);
  const double pi = std::acos(-1.0);
  const double scaled_fwhm = fwhm / slope;
  const double exp_factor = 4.0 * ln2 / (scaled_fwhm * scaled_fwhm);
  const double scale_factor = 2.0 * std::sqrt(ln2 / pi) / scaled_fwhm;
  std::vector<cplx> F(padded, cplx(0.0, 0.0));
  F[0] = cplx(scale_factor, 0.0);
  const usize half = padded / 2;
  for (usize n = 1; n <= half; ++n) {
    const double g = scale_factor * std::exp(-static_cast<double>(n * n) * exp_factor);
    F[n] = cplx(g, 0.0);
    F[padded - n] = cplx(g, 0.0);
  }
  std::vector<cplx> Ff(F);
  fft_inplace(Ff, false);

  // Wiener deconvolution of the histogram
  std::vector<cplx> Uf(padded);
  for (usize n = 0; n < padded; ++n) {
    const cplx c = std::conj(Ff[n]);
    const cplx G = c / (c * Ff[n] + noise);
    Uf[n] = Vf[n] * G.real();
  }
  std::vector<cplx> U(Uf);
  fft_inplace(U, true);
  for (usize n = 0; n < padded; ++n) U[n] = cplx(std::max(U[n].real(), 0.0), 0.0);

  // E(u | v) = (blur of u * U) / (blur of U)
  std::vector<cplx> numer(padded);
  for (usize n = 0; n < padded; ++n) {
    const double u = minv + (static_cast<double>(n) - static_cast<double>(hoff)) * slope;
    numer[n] = cplx(u * U[n].real(), 0.0);
  }
  fft_inplace(numer, false);
  for (usize n = 0; n < padded; ++n) numer[n] *= Ff[n];
  fft_inplace(numer, true);
  std::vector<cplx> denom(U);
  fft_inplace(denom, false);
  for (usize n = 0; n < padded; ++n) denom[n] *= Ff[n];
  fft_inplace(denom, true);

  std::vector<double> E(static_cast<usize>(bins), 0.0);
  for (int b = 0; b < bins; ++b) {
    const usize n = hoff + static_cast<usize>(b);
    const double d = denom[n].real();
    E[static_cast<usize>(b)] = (d != 0.0) ? numer[n].real() / d : 0.0;
  }

  for (usize v : idx) {
    const double cidx = (vals[v] - minv) / slope;
    int b = static_cast<int>(std::floor(cidx));
    if (b < 0) b = 0;
    double s;
    if (b < bins - 1) {
      s = E[static_cast<usize>(b)] + (E[static_cast<usize>(b + 1)] - E[static_cast<usize>(b)]) *
          (cidx - static_cast<double>(b));
    } else {
      s = E[static_cast<usize>(bins - 1)];
    }
    out[v] = s;
  }
}

// Coefficient of variation of exp(old - new) over the included voxels.
static double convergence_measure(const double* oldf, const double* newf,
                                  const std::vector<usize>& idx)
{
  const usize N = idx.size();
  if (N < 2) return 0.0;
  double mu = 0.0;
  for (usize v : idx) mu += std::exp(oldf[v] - newf[v]);
  mu /= static_cast<double>(N);
  double ss = 0.0;
  for (usize v : idx) {
    const double d = std::exp(oldf[v] - newf[v]) - mu;
    ss += d * d;
  }
  const double sd = std::sqrt(ss / static_cast<double>(N - 1));
  return sd / mu;
}

// ---------------------------------------------------------------------------
// 6. Grid geometry: ANTs padding rule and ITK shrink sampling
//
// Every quantity is computed in double and checked against explicit limits
// before it is converted to an integer, so that no conversion is out of range
// (undefined behaviour) and no lattice allocation can run away:
//   * the (virtual) padded grid may not exceed INT_MAX voxels along any axis,
//     so that padded indices stay representable on every platform;
//   * the control-point lattice may not exceed 2^27 points at the finest
//     level (about 1 GB per lattice copy), all axes together.
// All index arithmetic uses a 64-bit integer type (`long` is 32-bit on
// Windows).
// ---------------------------------------------------------------------------
typedef std::int64_t i64;

static const double kMaxPaddedSize = 2147483647.0;    // per axis
static const double kMaxControlPoints = 134217728.0;  // 2^27, finest lattice

struct AxisGeom {
  int n = 0;          // original size
  double spacing = 1; // mm
  int spans = 1;      // spline spans at level 0
  i64 lower = 0;      // padding added before index 0
  i64 padded = 0;     // padded size
  i64 out = 0;        // shrunk size (of the padded grid)
  i64 offset = 0;     // padded index of shrunk index 0
  i64 jlo = 0, jhi = -1; // shrunk indices that fall inside the original image
  int m = 0;          // jhi - jlo + 1
};

// `axis` (1-based) only names the dimension in error messages.
static AxisGeom axis_geometry(int n, double spacing, double distance, int shrink, int axis)
{
  AxisGeom g;
  g.n = n;
  g.spacing = spacing;
  const double domain = static_cast<double>(n - 1) * spacing;
  // number of spans covering the field of view (ANTs: ceil(domain / distance));
  // the relative tolerance makes an exact multiple, such as a distance equal
  // to the field of view, give exactly that number of spans despite round-off
  const double ratio = domain / distance;
  double spans_d = std::ceil(ratio - 1e-9 * ratio);
  if (!(spans_d >= 1.0)) spans_d = 1.0;
  if (!(spans_d <= kMaxControlPoints)) {
    Rcpp::stop("`n4_bias_field`: `spline_distance` (%g) is too small along dimension %d: it implies %g spline spans across a field of view of %g units. `spline_distance` is a distance in the units of `vox2ras` (usually millimeters), not a number of control points.",
               distance, axis, spans_d, domain);
  }
  double extra_d = std::floor((spans_d * distance - domain) / spacing + 0.5);
  if (!(extra_d > 0.0)) extra_d = 0.0;
  const double padded_d = static_cast<double>(n) + extra_d;
  if (!(padded_d <= kMaxPaddedSize)) {
    Rcpp::stop("`n4_bias_field`: `spline_distance` (%g) is too large along dimension %d: padding the field of view of %g units to a whole number of spline spans would need %g voxels (at most %g are supported); use a smaller `spline_distance`.",
               distance, axis, domain, padded_d, kMaxPaddedSize);
  }
  g.spans = static_cast<int>(spans_d);
  const i64 extra = static_cast<i64>(extra_d);
  g.lower = extra / 2;
  g.padded = static_cast<i64>(n) + extra;
  g.out = g.padded / shrink;
  if (g.out < 1) g.out = 1;
  const i64 rem = (g.padded - 1) - static_cast<i64>(shrink) * (g.out - 1);
  g.offset = (rem + 1) / 2;
  // shrunk index j samples padded index j * shrink + offset; keep those that
  // land on an original voxel (padded index in [lower, lower + n - 1])
  i64 jlo = 0;
  if (g.lower > g.offset) jlo = (g.lower - g.offset + shrink - 1) / shrink;
  i64 jhi = (g.lower + n - 1 - g.offset);
  jhi = (jhi < 0) ? -1 : jhi / shrink;
  if (jhi > g.out - 1) jhi = g.out - 1;
  g.jlo = jlo;
  g.jhi = jhi;
  g.m = static_cast<int>(jhi - jlo + 1);
  return g;
}

// Parametric coordinate of a padded index for a lattice with `spans` spans.
static inline double param_coord(i64 padded_index, i64 padded_size, int spans)
{
  if (padded_size <= 1) return 0.0;
  return static_cast<double>(spans) * static_cast<double>(padded_index) /
         static_cast<double>(padded_size - 1);
}

} // namespace raven4

// ===========================================================================
// Rcpp entry points (internal; wrapped by R/imaging-intensity-n4.R)
// ===========================================================================

// Lower / upper truncation bounds following the TruncateIntensity operation of
// ANTsPy 0.6.3 (verified against the installed binary): the histogram covers
// every finite voxel over [min, max], except that a minimum of exactly 0 is
// raised to 1e-6, so an exact-zero background (and anything below 1e-6) is
// left out. NA when no voxel falls into that range.
// [[Rcpp::export]]
Rcpp::NumericVector n4_truncate_quantiles(const Rcpp::NumericVector& x, double lower,
                                          double upper, int bins)
{
  using namespace raven4;
  if (bins < 1) Rcpp::stop("`n4_truncate_quantiles`: `bins` must be >= 1.");
  if (!(lower >= 0.0 && lower <= 1.0 && upper >= 0.0 && upper <= 1.0 && lower <= upper)) {
    Rcpp::stop("`n4_truncate_quantiles`: quantiles must satisfy 0 <= lower <= upper <= 1.");
  }
  const R_xlen_t n = x.size();
  double minv = std::numeric_limits<double>::infinity();
  double maxv = -std::numeric_limits<double>::infinity();
  for (R_xlen_t i = 0; i < n; ++i) {
    const double v = x[i];
    if (!std::isfinite(v)) continue;
    if (v < minv) minv = v;
    if (v > maxv) maxv = v;
  }
  Rcpp::NumericVector out(2);
  out[0] = NA_REAL; out[1] = NA_REAL;
  if (!(maxv >= minv)) return out;                  // no finite voxel
  const double lowb = (minv == 0.0) ? 1e-6 : minv;
  double count = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double v = x[i];
    if (std::isfinite(v) && v >= lowb) count += 1.0;
  }
  if (count == 0.0) return out;                     // e.g. an all-zero image
  if (!(maxv > lowb)) { out[0] = lowb; out[1] = maxv; return out; }
  std::vector<double> H(static_cast<usize>(bins), 0.0);
  const double interval = (maxv - lowb) / static_cast<double>(bins);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double v = x[i];
    if (!std::isfinite(v) || v < lowb) continue;
    int b = static_cast<int>(std::floor((v - lowb) / interval));
    if (b < 0) b = 0;
    if (b > bins - 1) b = bins - 1;
    H[static_cast<usize>(b)] += 1.0;
  }
  out[0] = histogram_quantile(H, count, lowb, maxv, lower);
  out[1] = histogram_quantile(H, count, lowb, maxv, upper);
  return out;
}

// Binary erosion (erode = TRUE) or dilation with the ITK ball of `radius`.
// [[Rcpp::export]]
Rcpp::LogicalVector n4_mask_morph(const Rcpp::LogicalVector& mask, const Rcpp::IntegerVector& dims,
                                  int radius, bool erode)
{
  using namespace raven4;
  check_dims(dims, mask.size(), "`n4_mask_morph`");
  if (radius < 0) Rcpp::stop("`n4_mask_morph`: `radius` must be >= 0.");
  std::vector<unsigned char> m = logical_to_bytes(mask);
  return bytes_to_logical(morph(m, dims[0], dims[1], dims[2], radius, erode));
}

// get_mask clean-up: erode(radius) -> [largest component] -> dilate(radius)
// -> fill holes. radius = 0 skips the erosion / dilation.
// [[Rcpp::export]]
Rcpp::LogicalVector n4_mask_cleanup(const Rcpp::LogicalVector& mask, const Rcpp::IntegerVector& dims,
                                    int radius, bool keep_largest, int min_size)
{
  using namespace raven4;
  check_dims(dims, mask.size(), "`n4_mask_cleanup`");
  if (radius < 0) Rcpp::stop("`n4_mask_cleanup`: `radius` must be >= 0.");
  if (min_size < 0) min_size = 0;
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  std::vector<unsigned char> m = logical_to_bytes(mask);
  if (radius > 0) m = morph(m, nx, ny, nz, radius, true);
  if (keep_largest) m = largest_component(m, nx, ny, nz, static_cast<usize>(min_size));
  if (radius > 0) m = morph(m, nx, ny, nz, radius, false);
  m = fill_holes(m, nx, ny, nz);
  return bytes_to_logical(m);
}

// Evaluate a control-point lattice on a regular grid covering its whole
// parametric domain (testing helper).
// [[Rcpp::export]]
Rcpp::NumericVector n4_bspline_evaluate(const Rcpp::NumericVector& lattice,
                                        const Rcpp::IntegerVector& ldims, int order,
                                        const Rcpp::IntegerVector& grid)
{
  using namespace raven4;
  if (order < 1 || order > 3) Rcpp::stop("`n4_bspline_evaluate`: `order` must be 1, 2 or 3.");
  check_dims(ldims, lattice.size(), "`n4_bspline_evaluate`");
  if (grid.size() != 3) Rcpp::stop("`n4_bspline_evaluate`: `grid` must have length 3.");
  Lattice lat;
  lat.alloc(ldims[0], ldims[1], ldims[2]);
  for (R_xlen_t i = 0; i < lattice.size(); ++i) lat.v[static_cast<usize>(i)] = lattice[i];
  AxisTable t[3];
  for (int d = 0; d < 3; ++d) {
    const int spans = ldims[d] - order;
    if (spans < 1) Rcpp::stop("`n4_bspline_evaluate`: lattice too small for this order.");
    if (grid[d] < 1) Rcpp::stop("`n4_bspline_evaluate`: invalid grid.");
    std::vector<double> u(static_cast<usize>(grid[d]));
    for (int i = 0; i < grid[d]; ++i) u[static_cast<usize>(i)] = param_coord(i, grid[d], spans);
    t[d] = make_axis_table(order, spans, u);
  }
  const R_xlen_t n = static_cast<R_xlen_t>(grid[0]) * grid[1] * grid[2];
  Rcpp::NumericVector out(n);
  reconstruct(lat, t[0], t[1], t[2], REAL(out));
  return out;
}

// Refine (double) a control-point lattice (testing helper).
// [[Rcpp::export]]
Rcpp::NumericVector n4_bspline_refine(const Rcpp::NumericVector& lattice,
                                      const Rcpp::IntegerVector& ldims, int order)
{
  using namespace raven4;
  if (order < 1 || order > 3) Rcpp::stop("`n4_bspline_refine`: `order` must be 1, 2 or 3.");
  check_dims(ldims, lattice.size(), "`n4_bspline_refine`");
  Lattice lat;
  lat.alloc(ldims[0], ldims[1], ldims[2]);
  for (R_xlen_t i = 0; i < lattice.size(); ++i) lat.v[static_cast<usize>(i)] = lattice[i];
  Lattice ref = refine_lattice(lat, order);
  Rcpp::NumericVector out(static_cast<R_xlen_t>(ref.size()));
  for (usize i = 0; i < ref.size(); ++i) out[static_cast<R_xlen_t>(i)] = ref.v[i];
  out.attr("dim") = Rcpp::IntegerVector::create(ref.nx, ref.ny, ref.nz);
  return out;
}

// The N4 core. Returns the log bias field on the full-resolution grid plus
// diagnostics. `mask` selects the voxels driving the fit (TRUE = include),
// `weight` is either empty or a per-voxel confidence (voxels with weight <= 0
// are excluded, the others weigh the B-spline fit).
// [[Rcpp::export]]
Rcpp::List n4_bias_field(const Rcpp::NumericVector& x, const Rcpp::IntegerVector& dims,
                         const Rcpp::LogicalVector& mask, const Rcpp::NumericVector& weight,
                         const Rcpp::NumericVector& spacing, int shrink,
                         const Rcpp::IntegerVector& iterations, double tol,
                         const Rcpp::NumericVector& spline_distance, int order, int bins,
                         double fwhm, double noise, bool verbose)
{
  using namespace raven4;
  check_dims(dims, x.size(), "`n4_bias_field`");
  if (mask.size() != x.size()) Rcpp::stop("`n4_bias_field`: `mask` must match the volume.");
  const bool use_weight = weight.size() > 0;
  if (use_weight && weight.size() != x.size()) Rcpp::stop("`n4_bias_field`: `weight` must match the volume.");
  if (spacing.size() != 3 || spline_distance.size() != 3) Rcpp::stop("`n4_bias_field`: `spacing` and `spline_distance` must have length 3.");
  if (shrink < 1) Rcpp::stop("`n4_bias_field`: `shrink` must be >= 1.");
  if (order < 1 || order > 3) Rcpp::stop("`n4_bias_field`: `order` must be 1, 2 or 3.");
  if (bins < 2) Rcpp::stop("`n4_bias_field`: `bins` must be >= 2.");
  if (!(fwhm > 0.0) || !(noise > 0.0)) Rcpp::stop("`n4_bias_field`: `fwhm` and `noise` must be > 0.");
  if (iterations.size() < 1) Rcpp::stop("`n4_bias_field`: `iterations` must have at least one level.");
  if (iterations.size() > 1024) Rcpp::stop("`n4_bias_field`: too many fitting levels (`iterations` has %d entries).", static_cast<int>(std::min<R_xlen_t>(iterations.size(), 2147483647)));
  for (R_xlen_t i = 0; i < iterations.size(); ++i) {
    if (iterations[i] == NA_INTEGER || iterations[i] < 0) Rcpp::stop("`n4_bias_field`: invalid `iterations`.");
  }
  for (int d = 0; d < 3; ++d) {
    if (dims[d] < 2) Rcpp::stop("`n4_bias_field`: every dimension must be at least 2.");
    if (!(spacing[d] > 0.0) || !std::isfinite(spacing[d])) Rcpp::stop("`n4_bias_field`: invalid `spacing`.");
    if (!(spline_distance[d] > 0.0) || !std::isfinite(spline_distance[d])) Rcpp::stop("`n4_bias_field`: invalid `spline_distance`.");
  }
  const int n_levels = static_cast<int>(iterations.size());

  // ---- geometry ----------------------------------------------------------
  AxisGeom geo[3];
  for (int d = 0; d < 3; ++d) {
    geo[d] = axis_geometry(dims[d], spacing[d], spline_distance[d], shrink, d + 1);
    if (geo[d].m < 1) {
      Rcpp::stop("`n4_bias_field`: the shrink factor (%d) leaves no sample along dimension %d; use a smaller `shrink_factor`.",
                 shrink, d + 1);
    }
  }
  // the lattice doubles its spans at every level: bound the finest one before
  // anything is allocated (computed in double, so any overflow shows as inf)
  double finest[3], finest_total = 1.0;
  for (int d = 0; d < 3; ++d) {
    finest[d] = std::ldexp(static_cast<double>(geo[d].spans), n_levels - 1) + static_cast<double>(order);
    finest_total *= finest[d];
  }
  if (!(finest_total <= kMaxControlPoints)) {
    Rcpp::stop("`n4_bias_field`: the B-spline lattice at the finest level would have %g x %g x %g control points (at most %g in total are supported); `spline_distance` is in the units of `vox2ras` (usually millimeters), so increase it or use fewer levels (`length(iterations)`).",
               finest[0], finest[1], finest[2], kMaxControlPoints);
  }
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const int mx = geo[0].m, my = geo[1].m, mz = geo[2].m;
  const usize Nb = static_cast<usize>(mx) * static_cast<usize>(my) * static_cast<usize>(mz);
  // a lattice much denser than the sub-sampled grid it is fitted to cannot be
  // estimated meaningfully; it usually means the distance is in the wrong units
  if (finest_total > 8.0 * static_cast<double>(Nb)) {
    Rcpp::warning("`n4_bias_field`: the B-spline lattice at the finest level has %g control points for a fitting grid of only %d x %d x %d samples; `spline_distance` is in the units of `vox2ras` (usually millimeters) and is probably far smaller than intended.",
                  finest_total, mx, my, mz);
  }

  // ---- shrunk box: included voxels, log intensities, weights -------------
  std::vector<unsigned char> inc(Nb, 0);
  std::vector<double> log_input(Nb, 0.0), wt(Nb, 1.0);
  std::vector<usize> inc_idx;
  {
    const usize sy = static_cast<usize>(nx), sz = static_cast<usize>(nx) * static_cast<usize>(ny);
    for (int jz = 0; jz < mz; ++jz) {
      const i64 iz = (geo[2].jlo + jz) * shrink + geo[2].offset - geo[2].lower;
      for (int jy = 0; jy < my; ++jy) {
        const i64 iy = (geo[1].jlo + jy) * shrink + geo[1].offset - geo[1].lower;
        for (int jx = 0; jx < mx; ++jx) {
          const i64 ix = (geo[0].jlo + jx) * shrink + geo[0].offset - geo[0].lower;
          const usize v = static_cast<usize>(ix) + sy * static_cast<usize>(iy) + sz * static_cast<usize>(iz);
          const usize b = static_cast<usize>(jx) + static_cast<usize>(mx) *
                          (static_cast<usize>(jy) + static_cast<usize>(my) * static_cast<usize>(jz));
          if (mask[static_cast<R_xlen_t>(v)] != TRUE) continue;
          const double val = x[static_cast<R_xlen_t>(v)];
          if (!std::isfinite(val)) continue;
          double w = 1.0;
          if (use_weight) {
            w = weight[static_cast<R_xlen_t>(v)];
            if (!std::isfinite(w) || !(w > 0.0)) continue;
          }
          inc[b] = 1;
          // as in ITK's N4 filter (and so ANTs / ANTsPy): a positive voxel
          // enters the fit through its logarithm, a zero or negative one with
          // its raw value
          log_input[b] = (val > 0.0) ? std::log(val) : val;
          wt[b] = w;
          inc_idx.push_back(b);
        }
      }
    }
  }
  if (inc_idx.empty()) {
    Rcpp::stop("`n4_bias_field`: no finite voxel inside the mask survives the shrinking; check `mask` / `weight_mask` or lower `shrink_factor`.");
  }

  std::vector<double> log_bias(Nb, 0.0), new_log_bias(Nb, 0.0);
  std::vector<double> log_uncorr(log_input), log_sharp(Nb, 0.0), residual(Nb, 0.0);

  // ---- lattice at level 0 --------------------------------------------------
  Lattice phi;
  phi.alloc(geo[0].spans + order, geo[1].spans + order, geo[2].spans + order);

  std::vector<double> trace;
  Rcpp::IntegerVector iterations_run(n_levels);

  if (verbose) {
    Rcpp::Rcout << "[N4] grid " << nx << " x " << ny << " x " << nz
                << ", padded " << geo[0].padded << " x " << geo[1].padded << " x " << geo[2].padded
                << ", shrunk box " << mx << " x " << my << " x " << mz
                << ", " << inc_idx.size() << " voxels drive the fit" << std::endl;
    Rcpp::Rcout << "[N4] level-0 control points " << phi.nx << " x " << phi.ny << " x " << phi.nz
                << " (spline order " << order << ")" << std::endl;
  }

  for (int level = 0; level < n_levels; ++level) {
    // parametric tables of the shrunk box for this level's lattice
    AxisTable tb[3];
    for (int d = 0; d < 3; ++d) {
      // bounded by the finest-lattice check above (level < 28, spans << level < 2^27)
      const int spans = static_cast<int>(static_cast<i64>(geo[d].spans) << level);
      std::vector<double> u(static_cast<usize>(geo[d].m));
      for (int j = 0; j < geo[d].m; ++j) {
        const i64 pj = (geo[d].jlo + j) * shrink + geo[d].offset;
        u[static_cast<usize>(j)] = param_coord(pj, geo[d].padded, spans);
      }
      tb[d] = make_axis_table(order, spans, u);
    }
    const usize Ls = phi.size();

    // fixed chunking of the box along z (independent of the thread count)
    usize G = static_cast<usize>(mz);
    if (G > 64) G = 64;
    const usize mem_cap = static_cast<usize>(4194304) / std::max<usize>(Ls, 1);  // ~64 MB of partials
    if (G > mem_cap) G = (mem_cap < 1) ? 1 : mem_cap;
    std::vector<int> chunk_z0(G), chunk_z1(G);
    for (usize g = 0; g < G; ++g) {
      chunk_z0[g] = static_cast<int>((g * static_cast<usize>(mz)) / G);
      chunk_z1[g] = static_cast<int>(((g + 1) * static_cast<usize>(mz)) / G);
    }
    std::vector<double> omega_part(G * Ls), delta_part(G * Ls);

    const int max_iter = iterations[level];
    if (verbose) {
      Rcpp::Rcout << "[N4] level " << (level + 1) << "/" << n_levels << ": control points "
                  << phi.nx << " x " << phi.ny << " x " << phi.nz << ", up to " << max_iter
                  << " iterations" << std::endl;
    }
    int it = 0;
    for (; it < max_iter; ++it) {
      // 1. sharpen the current estimate of the uncorrected log image
      sharpen_histogram(log_uncorr.data(), inc_idx, bins, fwhm, noise, log_sharp.data());
      // 2. residual log bias at the included voxels
      for (usize v : inc_idx) residual[v] = log_uncorr[v] - log_sharp[v];
      // 3. B-spline approximation of the residual; add it to the lattice
      std::fill(omega_part.begin(), omega_part.end(), 0.0);
      std::fill(delta_part.begin(), delta_part.end(), 0.0);
      {
        FitChunkWorker fw(inc.data(), residual.data(), wt.data(), mx, my,
                          &tb[0], &tb[1], &tb[2], phi.nx, phi.ny, Ls,
                          &chunk_z0, &chunk_z1, omega_part.data(), delta_part.data());
        TinyParallel::parallelFor(0, G, fw, 1);
        FitReduceWorker rw(omega_part.data(), delta_part.data(), G, Ls, phi.v.data());
        TinyParallel::parallelFor(0, Ls, rw, 2048);
      }
      // 4. smooth field on the box
      reconstruct(phi, tb[0], tb[1], tb[2], new_log_bias.data());
      // 5. convergence: CV of the ratio between successive bias estimates
      const double cv = convergence_measure(log_bias.data(), new_log_bias.data(), inc_idx);
      log_bias.swap(new_log_bias);
      for (usize v : inc_idx) log_uncorr[v] = log_input[v] - log_bias[v];
      trace.push_back(cv);
      if (verbose) {
        Rcpp::Rcout << "  iteration " << std::setw(3) << (it + 1) << " of " << max_iter
                    << ": convergence " << std::scientific << std::setprecision(4) << cv
                    << std::defaultfloat << " (threshold " << tol << ")" << std::endl;
      }
      Rcpp::checkUserInterrupt();
      if (cv <= tol) { ++it; break; }
    }
    iterations_run[level] = it;

    if (level < n_levels - 1) phi = refine_lattice(phi, order);
  }

  // ---- full-resolution reconstruction ----------------------------------------
  AxisTable tf[3];
  const int final_spans[3] = {phi.nx - order, phi.ny - order, phi.nz - order};
  for (int d = 0; d < 3; ++d) {
    std::vector<double> u(static_cast<usize>(dims[d]));
    for (int i = 0; i < dims[d]; ++i) u[static_cast<usize>(i)] = param_coord(i + geo[d].lower, geo[d].padded, final_spans[d]);
    tf[d] = make_axis_table(order, final_spans[d], u);
  }
  Rcpp::NumericVector out(x.size());
  reconstruct(phi, tf[0], tf[1], tf[2], REAL(out));

  return Rcpp::List::create(
    Rcpp::Named("log_bias") = out,
    Rcpp::Named("convergence") = Rcpp::wrap(trace),
    Rcpp::Named("iterations") = iterations_run,
    Rcpp::Named("control_points") = Rcpp::IntegerVector::create(phi.nx, phi.ny, phi.nz),
    Rcpp::Named("padded_dim") = Rcpp::IntegerVector::create(static_cast<int>(geo[0].padded), static_cast<int>(geo[1].padded), static_cast<int>(geo[2].padded)),
    Rcpp::Named("box_dim") = Rcpp::IntegerVector::create(mx, my, mz),
    Rcpp::Named("n_included") = static_cast<double>(inc_idx.size()));
}
