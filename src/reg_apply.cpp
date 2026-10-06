// Entry points for apply_transform3d_volume / apply_transform3d_points and the
// displacement-field inverse (see reg_apply.h for the shared machinery).
//
// A transform chain arrives from R as a list of stages, each a list with
// `type = "affine"` (plus `matrix`, 4x4 RAS->RAS) or `type = "field"` (plus
// `field`, an (nx, ny, nz, 3) RAS displacement array, `dim` and `vox2ras`).
// The chain is converted once, on the main thread, into borrowed views; the
// workers then stream every output voxel / point through it independently, so
// no composite field is ever materialized and the result does not depend on
// the thread count.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <string>
#include <cmath>
#include "reg_apply.h"

using namespace ravereg_apply;

namespace {

Matrix4d toMatrix4(const Rcpp::NumericMatrix& m) {
  if (m.nrow() != 4 || m.ncol() != 4) {
    Rcpp::stop("C++ `reg_apply`: expected a 4x4 matrix.");
  }
  Matrix4d M;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) M(i, j) = m(i, j);
  return M;
}

// Converted chain: the stage views plus the R vectors that own the field data.
struct ChainData {
  std::vector<Rcpp::NumericVector> keep;
  std::vector<Stage> stages;
};

void buildChain(const Rcpp::List& chain, ChainData& out) {
  const int n = chain.size();
  out.keep.reserve(n);
  out.stages.reserve(n);
  for (int s = 0; s < n; ++s) {
    Rcpp::List el = chain[s];
    const std::string type = Rcpp::as<std::string>(el["type"]);
    Stage st;
    if (type == "affine") {
      st.isField = false;
      st.M = toMatrix4(Rcpp::as<Rcpp::NumericMatrix>(el["matrix"]));
    } else if (type == "field") {
      st.isField = true;
      Rcpp::NumericVector f = Rcpp::as<Rcpp::NumericVector>(el["field"]);
      Rcpp::IntegerVector d = Rcpp::as<Rcpp::IntegerVector>(el["dim"]);
      if (d.size() < 3 || d[0] < 1 || d[1] < 1 || d[2] < 1) {
        Rcpp::stop("C++ `reg_apply`: a field stage needs three positive dimensions.");
      }
      st.nx = d[0]; st.ny = d[1]; st.nz = d[2];
      const std::size_t N = static_cast<std::size_t>(st.nx) *
                            static_cast<std::size_t>(st.ny) *
                            static_cast<std::size_t>(st.nz);
      if (static_cast<std::size_t>(f.size()) != 3 * N) {
        Rcpp::stop("C++ `reg_apply`: field stage %d has %d values, expected %d (nx * ny * nz * 3).",
                   s + 1, static_cast<int>(f.size()), static_cast<int>(3 * N));
      }
      out.keep.push_back(f);              // keeps the R object alive
      const double* base = f.begin();     // points into the R object itself
      st.ux = base; st.uy = base + N; st.uz = base + 2 * N;
      st.ras2vox = toMatrix4(Rcpp::as<Rcpp::NumericMatrix>(el["vox2ras"])).inverse();
    } else {
      Rcpp::stop("C++ `reg_apply`: unknown stage type '%s'.", type.c_str());
    }
    out.stages.push_back(st);
  }
}

// Resample every reference voxel (and every frame) through the chain.
struct VolumeWorker : public TinyParallel::Worker {
  const Stage* stages; int nStages;
  Matrix4d refV2R, volR2V;
  int rnx; std::size_t rsz, rN;              // reference grid: nx, nx*ny, nx*ny*nz
  const double* vol; int vnx, vny, vnz; std::size_t vN;
  int nFrames;
  int interp;                                // 0 = nearest, 1 = trilinear
  double naFill;
  double* out;

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t o = begin; o < end; ++o) {
      const std::size_t k = o / rsz;
      const std::size_t r = o - k * rsz;
      const std::size_t j = r / static_cast<std::size_t>(rnx);
      const std::size_t i = r - j * static_cast<std::size_t>(rnx);
      Vector3d p = (refV2R * Vector4d(static_cast<double>(i), static_cast<double>(j),
                                      static_cast<double>(k), 1.0)).head<3>();
      bool ok = mapPoint(stages, nStages, p);
      double cx = 0.0, cy = 0.0, cz = 0.0;
      if (ok) {
        const Vector4d c = volR2V * Vector4d(p[0], p[1], p[2], 1.0);
        cx = c[0]; cy = c[1]; cz = c[2];
        ok = std::isfinite(cx) && std::isfinite(cy) && std::isfinite(cz);
      }
      if (!ok) {
        for (int f = 0; f < nFrames; ++f) out[o + rN * f] = naFill;
        continue;
      }
      if (interp == 0) {
        // nearest neighbour: round, then bounds-check the index (resample3D's rule);
        // compared as doubles so an enormous coordinate never overflows a cast
        const double rx = std::nearbyint(cx), ry = std::nearbyint(cy), rz = std::nearbyint(cz);
        if (rx < 0.0 || ry < 0.0 || rz < 0.0 ||
            rx >= static_cast<double>(vnx) || ry >= static_cast<double>(vny) ||
            rz >= static_cast<double>(vnz)) {
          for (int f = 0; f < nFrames; ++f) out[o + rN * f] = naFill;
          continue;
        }
        const std::size_t idx = static_cast<std::size_t>(rx) +
          static_cast<std::size_t>(vnx) * (static_cast<std::size_t>(ry) +
          static_cast<std::size_t>(vny) * static_cast<std::size_t>(rz));
        for (int f = 0; f < nFrames; ++f) out[o + rN * f] = vol[idx + vN * f];
      } else {
        for (int f = 0; f < nFrames; ++f) {
          double v;
          const bool in = ravereg::trilinearSample(vol + vN * f, vnx, vny, vnz, cx, cy, cz, v);
          out[o + rN * f] = in ? v : naFill;
        }
      }
    }
  }
};

// Map every point (row of a column-major N x 3 matrix) through the chain.
struct PointsWorker : public TinyParallel::Worker {
  const Stage* stages; int nStages;
  const double* in; double* out; std::size_t n;
  double naValue;

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t r = begin; r < end; ++r) {
      Vector3d p(in[r], in[r + n], in[r + 2 * n]);
      const bool ok = std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) &&
                      mapPoint(stages, nStages, p);
      if (ok) {
        out[r] = p[0]; out[r + n] = p[1]; out[r + 2 * n] = p[2];
      } else {
        out[r] = naValue; out[r + n] = naValue; out[r + 2 * n] = naValue;
      }
    }
  }
};

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::NumericVector apply_transform3d_volume_cpp(const Rcpp::NumericVector& volume,
                                                 const Rcpp::IntegerVector& volumeDim,
                                                 const Rcpp::NumericMatrix& volumeVox2Ras,
                                                 const Rcpp::IntegerVector& referenceDim,
                                                 const Rcpp::NumericMatrix& referenceVox2Ras,
                                                 const Rcpp::List& chain,
                                                 const int interpolation,
                                                 const double naFill) {
  if (volumeDim.size() < 3 || referenceDim.size() < 3) {
    Rcpp::stop("C++ `apply_transform3d_volume_cpp`: dimensions must have length 3 (or 4 for frames).");
  }
  const int vnx = volumeDim[0], vny = volumeDim[1], vnz = volumeDim[2];
  const int nFrames = (volumeDim.size() >= 4) ? volumeDim[3] : 1;
  const int rnx = referenceDim[0], rny = referenceDim[1], rnz = referenceDim[2];
  if (vnx < 1 || vny < 1 || vnz < 1 || nFrames < 1 || rnx < 1 || rny < 1 || rnz < 1) {
    Rcpp::stop("C++ `apply_transform3d_volume_cpp`: dimensions must be positive.");
  }
  const std::size_t vN = static_cast<std::size_t>(vnx) * static_cast<std::size_t>(vny) *
                         static_cast<std::size_t>(vnz);
  if (static_cast<std::size_t>(volume.size()) != vN * static_cast<std::size_t>(nFrames)) {
    Rcpp::stop("C++ `apply_transform3d_volume_cpp`: `volume` length does not match its dimensions.");
  }
  const std::size_t rN = static_cast<std::size_t>(rnx) * static_cast<std::size_t>(rny) *
                         static_cast<std::size_t>(rnz);

  ChainData cd;
  buildChain(chain, cd);

  Rcpp::NumericVector out(static_cast<R_xlen_t>(rN * static_cast<std::size_t>(nFrames)));

  VolumeWorker w;
  w.stages = cd.stages.data(); w.nStages = static_cast<int>(cd.stages.size());
  w.refV2R = toMatrix4(referenceVox2Ras);
  w.volR2V = toMatrix4(volumeVox2Ras).inverse();
  w.rnx = rnx; w.rsz = static_cast<std::size_t>(rnx) * static_cast<std::size_t>(rny); w.rN = rN;
  w.vol = volume.begin(); w.vnx = vnx; w.vny = vny; w.vnz = vnz; w.vN = vN;
  w.nFrames = nFrames;
  w.interp = (interpolation == 0) ? 0 : 1;
  w.naFill = naFill;
  w.out = out.begin();
  TinyParallel::parallelFor(0, rN, w, 256);

  if (volumeDim.size() >= 4) {
    out.attr("dim") = Rcpp::IntegerVector::create(rnx, rny, rnz, nFrames);
  } else {
    out.attr("dim") = Rcpp::IntegerVector::create(rnx, rny, rnz);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix apply_transform3d_points_cpp(const Rcpp::NumericMatrix& points,
                                                 const Rcpp::List& chain) {
  if (points.ncol() != 3) {
    Rcpp::stop("C++ `apply_transform3d_points_cpp`: `points` must have 3 columns.");
  }
  const std::size_t n = static_cast<std::size_t>(points.nrow());
  ChainData cd;
  buildChain(chain, cd);

  Rcpp::NumericMatrix out(points.nrow(), 3);
  if (n == 0) return out;

  PointsWorker w;
  w.stages = cd.stages.data(); w.nStages = static_cast<int>(cd.stages.size());
  w.in = points.begin(); w.out = out.begin(); w.n = n;
  w.naValue = NA_REAL;
  TinyParallel::parallelFor(0, n, w, 1024);
  return out;
}

// [[Rcpp::export]]
Rcpp::List invert_displacement_field_cpp(const Rcpp::NumericVector& field,
                                         const Rcpp::IntegerVector& dim,
                                         const Rcpp::NumericMatrix& vox2ras,
                                         const int maxIterations,
                                         const double meanTolerance,
                                         const double maxTolerance) {
  if (dim.size() < 3 || dim[0] < 1 || dim[1] < 1 || dim[2] < 1) {
    Rcpp::stop("C++ `invert_displacement_field_cpp`: `dim` must hold three positive integers.");
  }
  const int nx = dim[0], ny = dim[1], nz = dim[2];
  const std::size_t N = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                        static_cast<std::size_t>(nz);
  if (static_cast<std::size_t>(field.size()) != 3 * N) {
    Rcpp::stop("C++ `invert_displacement_field_cpp`: `field` must hold nx * ny * nz * 3 values.");
  }
  const double* u = field.begin();
  Rcpp::NumericVector out(static_cast<R_xlen_t>(3 * N));
  double* v = out.begin();

  InvertOptions opt;
  opt.maxIterations = maxIterations;
  opt.meanTolerance = meanTolerance;
  opt.maxTolerance = maxTolerance;
  const InvertResult res = invertDisplacementField<double>(
    u, u + N, u + 2 * N, nx, ny, nz, toMatrix4(vox2ras), v, v + N, v + 2 * N, opt);

  out.attr("dim") = Rcpp::IntegerVector::create(nx, ny, nz, 3);
  return Rcpp::List::create(
    Rcpp::Named("field") = out,
    Rcpp::Named("iterations") = res.iterations,
    Rcpp::Named("mean_residual") = res.meanResidual,
    Rcpp::Named("max_residual") = res.maxResidual);
}
