// Prior-guided Gaussian-mixture tissue segmentation with a mean-field Markov
// random field ("Atropos-lite"), the native counterpart of ANTs Atropos.
//
// Model (K classes, voxel i inside the mask with intensity y_i), following the
// default ("Socrates") posterior formulation of Atropos:
//
//   P_i(k) proportional to  s_k(i)^w * ( N(y_i | mu_k, sigma_k^2) * M_k(i) )^(1 - w)
//
//   s_k(i) = pi_k * p_k(i) / sum_c pi_c * p_c(i)             spatial prior
//   M_k(i) = exp(beta * S_k(i)) / sum_c exp(beta * S_c(i))   MRF term
//   S_k(i) = sum_j omega_ij * P_j(k)                         neighborhood sum
//
// where w is the prior weight, pi_k the mixture proportions, p_k the per-voxel
// normalized spatial prior, beta the MRF smoothing factor and omega_ij = 1 /
// (physical distance between voxels i and j) over the (2 r + 1)^3 - 1
// neighborhood. As in Atropos, s_k, M_k and the Gaussian density are each
// floored at 1e-10 before the powers are taken. Without priors the spatial
// prior is the same constant for every class and Atropos applies no prior
// weight, so P_i(k) is proportional to N * M_k and the proportions do not enter
// the posterior. With w = 0 the priors only initialize the model; with w = 1 the
// posterior is s_k(i) itself and the image is ignored.
//
// The MRF term is the mean-field form of the Potts prior of Atropos: the
// neighbors' posteriors P_j(k) replace the hard labels that Atropos plugs in
// (after an iterated-conditional-modes sweep). Each iteration is one
// synchronous E-step (all voxels read the previous posteriors; initially the
// one-hot initial labels) followed by an M-step that re-estimates mu, sigma
// (posterior-weighted mean and unbiased weighted variance, the ITK estimators
// Atropos uses, with a variance floor) and pi (mean posterior) from the new
// posteriors.
//
// The original plan formula  N(y_i | mu_k, sigma_k^2) * pi_k^(1 - w) * p_k(i)^w
// * exp(beta * S_k(i))  is kept as commented-out code right next to the active
// formula in GmmEStepWorker::operator() ("ALTERNATIVE FORMULA").
//
// Initialization: either deterministic k-means on the masked intensities
// (seeded at equally spaced intensities, classes ordered by ascending mean) or
// the argmax of the priors (classes follow the prior order; voxels where every
// prior is zero start unlabeled). The initial Gaussian parameters come from
// the hard initial labeling, weighted by the prior value in the prior case,
// and the initial proportions are the label fractions, as in Atropos.
//
// Convergence (the rule of Atropos): after every iteration the mean maximum
// posterior over the mask is compared with the previous iteration's value and
// the loop stops as soon as the signed change is below `tolerance` (the first
// iteration never stops), so with tolerance = 0 the first decrease stops it;
// the posteriors and labels of the iteration just completed are returned.
//
// Implementation notes:
//   * Computation is restricted to the mask: intensities, priors and the two
//     posterior buffers are stored compactly (O(K * N_mask)); a full-grid
//     int map (grid index -> compact index, -1 outside) serves the neighbor
//     look-ups. Nothing outside the mask is ever read.
//   * The E-step runs in parallel over fixed-size blocks of mask voxels. Every
//     block's partial sums (for the M-step and the mean maximum posterior) are
//     written to a block-indexed array and combined serially in block order,
//     so the result is bit-identical regardless of the thread count.
//   * Worker threads never touch the R API; interrupts are checked between
//     iterations in the serial driver.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstddef>
#include <iomanip>
#include "TinyParallel.h"

namespace {

typedef std::size_t usize;

// Probability floor applied to the spatial prior, the MRF term and the
// Gaussian density before the powers are taken (Atropos uses the same 1e-10).
const double kProbEps = 1e-10;

// Voxels per reduction block. Fixed, so the partial sums do not depend on how
// TinyParallel splits the block range among threads.
const usize kBlock = 2048;

// Maximum Lloyd iterations for the k-means initialization (ANTs uses 200).
const int kKMeansMaxIter = 200;

struct Neighbor {
  int dx, dy, dz;
  std::ptrdiff_t delta;   // flat-index offset on the full grid
  double weight;          // 1 / physical distance
};

inline void decode_index(usize o, usize nx, usize nxy, int& i, int& j, int& k) {
  const usize kk = o / nxy;
  const usize r = o - kk * nxy;
  const usize jj = r / nx;
  const usize ii = r - jj * nx;
  i = static_cast<int>(ii);
  j = static_cast<int>(jj);
  k = static_cast<int>(kk);
}

// Unbiased weighted variance from the weighted sums (ITK's weighted covariance
// estimator used by Atropos): sum w (y - mean)^2 / (sum w - sum w^2 / sum w).
// `ss` is the centered weighted sum of squares, `s0` the sum of the weights and
// `sq` the sum of the squared weights.
inline double weighted_variance(double ss, double s0, double sq, double floor) {
  const double denom = s0 - sq / s0;
  if (!(denom > 0.0) || !(ss > 0.0)) return floor;
  return std::max(ss / denom, floor);
}

// One synchronous E-step over blocks [b_begin, b_end) of mask voxels.
struct GmmEStepWorker : public TinyParallel::Worker {
  const double* y;                 // compact intensities (N)
  const usize* grid;               // compact -> full-grid flat index (N)
  const int* map;                  // full grid -> compact index, -1 outside
  int nx, ny, nz;
  usize nxy;
  const std::vector<Neighbor>* nbrs;
  const double* post_old;          // N x K, voxel-major (previous posteriors)
  const double* prior;             // N x K normalized priors, or nullptr (w == 0)
  const double* mu;                // K
  const double* log_norm;          // K: -0.5 * log(2 pi sigma^2)
  const double* half_inv_var;      // K: 0.5 / sigma^2
  const double* pi;                // K mixture proportions
  const double* log_pi;            // K: log(max(pi, eps)); alternative formula only
  int K;
  usize N;
  double beta;
  double w;                        // prior weight (0 without priors)
  double* post_new;                // N x K, voxel-major (output)
  double* block_sums;              // n_blocks x (4 K + 1)

  GmmEStepWorker()
    : y(nullptr), grid(nullptr), map(nullptr), nx(0), ny(0), nz(0), nxy(0),
      nbrs(nullptr), post_old(nullptr), prior(nullptr), mu(nullptr),
      log_norm(nullptr), half_inv_var(nullptr), pi(nullptr), log_pi(nullptr),
      K(0), N(0), beta(0.0), w(0.0), post_new(nullptr), block_sums(nullptr) {}

  void operator()(usize b_begin, usize b_end) override {
    const usize Ku = static_cast<usize>(K);
    const bool use_mrf = (beta > 0.0) && nbrs && !nbrs->empty();
    const double log_eps = std::log(kProbEps);
    const double one_minus_w = 1.0 - w;
    // mrf: neighborhood sums S_k (stay 0 without the MRF); log_mrf: log M_k
    std::vector<double> mrf(Ku, 0.0), log_mrf(Ku, 0.0), sp(Ku, 0.0), lp(Ku, 0.0);
    const usize n_sums = 4 * Ku + 1;

    for (usize b = b_begin; b < b_end; ++b) {
      const usize m0 = b * kBlock;
      const usize m1 = std::min(N, m0 + kBlock);
      double* bs = block_sums + b * n_sums;
      for (usize s = 0; s < n_sums; ++s) bs[s] = 0.0;

      for (usize m = m0; m < m1; ++m) {
        // Mean-field neighborhood sums S_k(i) = sum_j omega_ij P_j(k) and the
        // MRF term M_k = exp(beta S_k) / sum_c exp(beta S_c), floored, in logs
        if (use_mrf) {
          for (usize k = 0; k < Ku; ++k) mrf[k] = 0.0;
          int i, j, kz;
          const usize og = grid[m];
          decode_index(og, static_cast<usize>(nx), nxy, i, j, kz);
          const std::vector<Neighbor>& nb = *nbrs;
          for (usize t = 0; t < nb.size(); ++t) {
            const int ii = i + nb[t].dx;
            if (ii < 0 || ii >= nx) continue;
            const int jj = j + nb[t].dy;
            if (jj < 0 || jj >= ny) continue;
            const int kk = kz + nb[t].dz;
            if (kk < 0 || kk >= nz) continue;
            const usize oj = static_cast<usize>(static_cast<std::ptrdiff_t>(og) + nb[t].delta);
            const int mj = map[oj];
            if (mj < 0) continue;
            const double* pj = post_old + static_cast<usize>(mj) * Ku;
            const double wgt = nb[t].weight;
            for (usize k = 0; k < Ku; ++k) mrf[k] += wgt * pj[k];
          }
          double mx = -std::numeric_limits<double>::infinity();
          for (usize k = 0; k < Ku; ++k) {
            log_mrf[k] = beta * mrf[k];
            if (log_mrf[k] > mx) mx = log_mrf[k];
          }
          double s = 0.0;
          for (usize k = 0; k < Ku; ++k) s += std::exp(log_mrf[k] - mx);
          const double lse = mx + std::log(s);
          for (usize k = 0; k < Ku; ++k) log_mrf[k] = std::max(log_mrf[k] - lse, log_eps);
        }

        // Spatial prior s_k(i): the normalized priors re-weighted by the
        // mixture proportions; where every prior is zero it stays zero (and
        // the floor below makes it uniform), as in Atropos
        if (prior) {
          const double* pm = prior + m * Ku;
          double den = 0.0;
          for (usize k = 0; k < Ku; ++k) den += pi[k] * pm[k];
          for (usize k = 0; k < Ku; ++k) sp[k] = (den > 0.0) ? pm[k] * pi[k] / den : pm[k];
        }

        const double yv = y[m];
        double mx = -std::numeric_limits<double>::infinity();
        for (usize k = 0; k < Ku; ++k) {
          const double d = yv - mu[k];
          const double log_lik_raw = log_norm[k] - half_inv_var[k] * d * d;
          const double log_lik = std::max(log_lik_raw, log_eps);
          const double log_spatial = prior ? std::log(std::max(sp[k], kProbEps)) : log_eps;
          // Atropos (Socrates): s_k^w * (N(y | mu_k, sigma_k^2) * M_k)^(1 - w)
          double v = w * log_spatial + one_minus_w * (log_lik + log_mrf[k]);
          // ALTERNATIVE FORMULA (original plan): N(y | mu_k, sigma_k^2) * pi_k^(1 - w)
          // * p_k^w * exp(beta * S_k). To test it, comment out the assignment
          // above and uncomment the next line (no other change is needed):
          // double v = log_lik_raw + one_minus_w * log_pi[k] + w * std::log(std::max(prior ? prior[m * Ku + k] : 1.0, kProbEps)) + beta * mrf[k];
          lp[k] = v;
          if (v > mx) mx = v;
        }
        double s = 0.0;
        for (usize k = 0; k < Ku; ++k) {
          lp[k] = std::exp(lp[k] - mx);
          s += lp[k];
        }
        const double inv = 1.0 / s;
        double pmax = 0.0;
        double* pn = post_new + m * Ku;
        for (usize k = 0; k < Ku; ++k) {
          const double p = lp[k] * inv;
          pn[k] = p;
          if (p > pmax) pmax = p;
          const double d = yv - mu[k];
          bs[4 * k] += p;
          bs[4 * k + 1] += p * d;
          bs[4 * k + 2] += p * d * d;
          bs[4 * k + 3] += p * p;
        }
        bs[4 * Ku] += pmax;
      }
    }
  }
};

// Deterministic 1D k-means (Lloyd) with equally spaced seeds; returns the
// centroids sorted ascending. Empty clusters keep their centroid. The cluster
// sums accumulate intensities relative to `ymin`, so a large common offset
// (intensities of the form offset + signal) does not swamp the signal in the
// running sums.
std::vector<double> kmeans_1d(const std::vector<double>& y, int K,
                              double ymin, double ymax) {
  const usize Ku = static_cast<usize>(K);
  std::vector<double> c(Ku);
  for (usize k = 0; k < Ku; ++k) {
    c[k] = ymin + (ymax - ymin) * (static_cast<double>(k) + 0.5) / static_cast<double>(K);
  }
  std::vector<double> sum(Ku);
  std::vector<usize> cnt(Ku);
  for (int it = 0; it < kKMeansMaxIter; ++it) {
    std::fill(sum.begin(), sum.end(), 0.0);
    std::fill(cnt.begin(), cnt.end(), static_cast<usize>(0));
    for (usize m = 0; m < y.size(); ++m) {
      const double v = y[m];
      usize best = 0;
      double bd = std::abs(v - c[0]);
      for (usize k = 1; k < Ku; ++k) {
        const double d = std::abs(v - c[k]);
        if (d < bd) { bd = d; best = k; }
      }
      sum[best] += v - ymin;
      cnt[best] += 1;
    }
    bool changed = false;
    for (usize k = 0; k < Ku; ++k) {
      if (cnt[k] > 0) {
        const double nc = ymin + sum[k] / static_cast<double>(cnt[k]);
        if (nc != c[k]) changed = true;
        c[k] = nc;
      }
    }
    Rcpp::checkUserInterrupt();
    if (!changed) break;
  }
  std::sort(c.begin(), c.end());
  return c;
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List segment_gmm_mrf_cpp(const Rcpp::NumericVector& volume,
                               const Rcpp::IntegerVector& dims,
                               const Rcpp::LogicalVector& mask,
                               const Rcpp::List& priors,
                               int n_classes,
                               double prior_weight,
                               double mrf_beta,
                               const Rcpp::IntegerVector& mrf_radius,
                               const Rcpp::NumericMatrix& direction,
                               int iterations,
                               double tolerance,
                               bool verbose) {
  if (dims.size() != 3) Rcpp::stop("`dims` must have length 3.");
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  if (nx <= 0 || ny <= 0 || nz <= 0) Rcpp::stop("`dims` must be positive.");
  const usize nxy = static_cast<usize>(nx) * static_cast<usize>(ny);
  const usize n_full = nxy * static_cast<usize>(nz);
  if (static_cast<usize>(volume.size()) != n_full || static_cast<usize>(mask.size()) != n_full) {
    Rcpp::stop("`volume` and `mask` must have prod(dims) elements.");
  }
  if (n_classes < 2) Rcpp::stop("`n_classes` must be at least 2.");
  const int K = n_classes;
  const usize Ku = static_cast<usize>(K);
  const bool use_priors = priors.size() > 0;
  if (use_priors && priors.size() != K) {
    Rcpp::stop("The number of priors must equal the number of classes.");
  }
  if (mrf_radius.size() != 3) Rcpp::stop("`mrf_radius` must have length 3.");
  if (direction.nrow() != 3 || direction.ncol() != 3) {
    Rcpp::stop("`direction` must be a 3 x 3 matrix.");
  }
  if (iterations < 1) Rcpp::stop("`iterations` must be at least 1.");
  if (prior_weight < 0.0 || prior_weight > 1.0) Rcpp::stop("`prior_weight` must be in [0, 1].");
  if (mrf_beta < 0.0) Rcpp::stop("`mrf_beta` must be non-negative.");
  // Without prior images Atropos applies no prior weight (the spatial prior is
  // the same constant for every class)
  const double w = use_priors ? prior_weight : 0.0;

  // ---- compact mask layout -------------------------------------------------
  std::vector<int> map(n_full, -1);
  usize N = 0;
  for (usize o = 0; o < n_full; ++o) {
    if (mask[o] == TRUE) ++N;
  }
  if (N == 0) Rcpp::stop("The mask contains no voxels.");
  if (N < Ku) Rcpp::stop("The mask contains fewer voxels than classes.");
  if (N > static_cast<usize>(std::numeric_limits<int>::max())) {
    Rcpp::stop("The mask contains too many voxels.");
  }
  std::vector<usize> grid(N);
  std::vector<double> y(N);
  {
    usize m = 0;
    for (usize o = 0; o < n_full; ++o) {
      if (mask[o] != TRUE) continue;
      const double v = volume[o];
      if (!std::isfinite(v)) Rcpp::stop("Intensities inside the mask must be finite.");
      map[o] = static_cast<int>(m);
      grid[m] = o;
      y[m] = v;
      ++m;
    }
  }

  // Global statistics of the masked intensities (variance floor, k-means
  // seeds). The mean is accumulated relative to the minimum so that a large
  // common offset does not swamp the signal in the running sum.
  double ymin = y[0], ymax = y[0];
  for (usize m = 0; m < N; ++m) {
    if (y[m] < ymin) ymin = y[m];
    if (y[m] > ymax) ymax = y[m];
  }
  double ysum = 0.0;
  for (usize m = 0; m < N; ++m) ysum += y[m] - ymin;
  const double gmean = ymin + ysum / static_cast<double>(N);
  double gvar = 0.0;
  for (usize m = 0; m < N; ++m) {
    const double d = y[m] - gmean;
    gvar += d * d;
  }
  gvar /= static_cast<double>(N);
  if (!std::isfinite(gvar)) {
    // squared deviations overflow (intensity range beyond about 1e154): every
    // downstream density would be NaN
    Rcpp::stop("Intensities inside the mask span too large a range (%g to %g) for their variance to be computed; rescale the volume first.",
               ymin, ymax);
  }
  if (!(gvar > 0.0) || !(ymax > ymin)) {
    Rcpp::stop("Intensities inside the mask are constant; nothing to segment.");
  }
  const double var_floor = 1e-6 * gvar;

  // ---- initialization ------------------------------------------------------
  std::vector<int> label(N, 0);          // initial hard labels (0-based; -1 = unlabeled)
  std::vector<double> init_w(N, 1.0);    // weight of each voxel in its class
  std::vector<double> prior;             // N x K normalized priors, only when w > 0
  // soft fallback statistics (prior-weighted over every voxel) for classes
  // that are never the argmax of the priors
  std::vector<double> soft_s0(Ku, 0.0), soft_s1(Ku, 0.0), soft_s2(Ku, 0.0), soft_sq(Ku, 0.0);

  if (use_priors) {
    // A list element that is not already a double vector (integer or logical
    // priors) is coerced to a fresh double copy that only its NumericVector
    // handle protects from the garbage collector, so every handle must stay
    // alive for as long as its data pointer is dereferenced below.
    std::vector<Rcpp::NumericVector> prior_vec(Ku);
    std::vector<const double*> pp(Ku);
    for (usize k = 0; k < Ku; ++k) {
      Rcpp::NumericVector p = priors[k];
      if (static_cast<usize>(p.size()) != n_full) {
        Rcpp::stop("Prior %d must have prod(dims) elements.", static_cast<int>(k + 1));
      }
      prior_vec[k] = p;
      pp[k] = REAL(prior_vec[k]);
    }
    if (w > 0.0) prior.assign(N * Ku, 0.0);
    std::vector<double> pv(Ku);
    usize n_informative = 0;   // mask voxels where at least one prior is positive
    for (usize m = 0; m < N; ++m) {
      const usize o = grid[m];
      double s = 0.0;
      for (usize k = 0; k < Ku; ++k) {
        const double v = pp[k][o];
        if (!(v >= 0.0) || !std::isfinite(v)) {
          Rcpp::stop("Prior %d has a negative or non-finite value inside the mask.",
                     static_cast<int>(k + 1));
        }
        pv[k] = v;
        s += v;
      }
      if (s > 0.0) {
        ++n_informative;
        usize best = 0;
        for (usize k = 0; k < Ku; ++k) {
          pv[k] /= s;
          if (pv[k] > pv[best]) best = k;
        }
        label[m] = static_cast<int>(best);
        init_w[m] = pv[best];
      } else {
        // Every prior is zero: the voxel starts unlabeled (label 0 in Atropos)
        // and its spatial prior stays zero (uniform after the floor)
        for (usize k = 0; k < Ku; ++k) pv[k] = 0.0;
        label[m] = -1;
        init_w[m] = 0.0;
      }
      const double d = y[m] - gmean;
      for (usize k = 0; k < Ku; ++k) {
        soft_s0[k] += pv[k];
        soft_s1[k] += pv[k] * d;
        soft_s2[k] += pv[k] * d * d;
        soft_sq[k] += pv[k] * pv[k];
        if (w > 0.0) prior[m * Ku + k] = pv[k];
      }
    }
    if (n_informative == 0) {
      // Every class would start from the same (global) statistics and the EM
      // symmetry would never be broken: a one-class result with duplicated
      // classes is never what the caller wants.
      Rcpp::stop("The priors are zero at every voxel inside the mask and carry no information; check that they are probability maps on the grid of `volume` (not integer-truncated) and that they cover the mask.");
    }
  } else {
    const std::vector<double> c = kmeans_1d(y, K, ymin, ymax);
    for (usize m = 0; m < N; ++m) {
      const double v = y[m];
      usize best = 0;
      double bd = std::abs(v - c[0]);
      for (usize k = 1; k < Ku; ++k) {
        const double d = std::abs(v - c[k]);
        if (d < bd) { bd = d; best = k; }
      }
      label[m] = static_cast<int>(best);
    }
  }

  // Initial Gaussian parameters and proportions from the hard labeling
  std::vector<double> mu(Ku, 0.0), var(Ku, 0.0), pi(Ku, 0.0);
  {
    std::vector<double> s0(Ku, 0.0), s1(Ku, 0.0), s2(Ku, 0.0), sq(Ku, 0.0);
    std::vector<usize> cnt(Ku, 0);
    usize n_labeled = 0;
    for (usize m = 0; m < N; ++m) {
      if (label[m] < 0) continue;
      const usize k = static_cast<usize>(label[m]);
      const double wt = init_w[m];
      const double d = y[m] - gmean;
      s0[k] += wt;
      s1[k] += wt * d;
      s2[k] += wt * d * d;
      sq[k] += wt * wt;
      cnt[k] += 1;
      ++n_labeled;
    }
    for (usize k = 0; k < Ku; ++k) {
      double a0 = s0[k], a1 = s1[k], a2 = s2[k], aq = sq[k];
      if (!(a0 > 0.0) && use_priors) {
        a0 = soft_s0[k]; a1 = soft_s1[k]; a2 = soft_s2[k]; aq = soft_sq[k];
      }
      if (!(a0 > 0.0)) {
        if (use_priors) {
          Rcpp::stop("Class %d received no voxels during initialization: prior %d is zero at every voxel inside the mask; check the priors and the mask.",
                     static_cast<int>(k + 1), static_cast<int>(k + 1));
        }
        // The k-means seeds are spread evenly over [min, max] of the masked
        // intensities, so a few extreme values (hot voxels, untruncated CT)
        // leave the middle seeds without any voxel.
        Rcpp::stop("Class %d received no voxels during the k-means initialization: the intensity range inside the mask (%g to %g) is dominated by extreme values. Truncate or clip the intensities (for example with the intensity truncation of `bias_correction_n4`), tighten the mask, or reduce `n_classes`.",
                   static_cast<int>(k + 1), ymin, ymax);
      }
      const double e1 = a1 / a0;
      mu[k] = gmean + e1;
      var[k] = weighted_variance(a2 - e1 * a1, a0, aq, var_floor);
      pi[k] = static_cast<double>(cnt[k]) / static_cast<double>(n_labeled);
    }
  }

  // ---- MRF neighborhood ----------------------------------------------------
  std::vector<Neighbor> nbrs;
  if (mrf_beta > 0.0) {
    if (mrf_radius[0] < 0 || mrf_radius[1] < 0 || mrf_radius[2] < 0) {
      Rcpp::stop("`mrf_radius` must be non-negative.");
    }
    // An offset of |d| >= n along an axis can never land on a voxel of the
    // grid, so clamping the radius to n - 1 cannot change any result; it only
    // keeps the neighbor list (and the uninterruptible E-step) bounded.
    const int rx = std::min(mrf_radius[0], nx - 1);
    const int ry = std::min(mrf_radius[1], ny - 1);
    const int rz = std::min(mrf_radius[2], nz - 1);
    for (int dz = -rz; dz <= rz; ++dz) {
      for (int dy = -ry; dy <= ry; ++dy) {
        for (int dx = -rx; dx <= rx; ++dx) {
          if (dx == 0 && dy == 0 && dz == 0) continue;
          double dist2 = 0.0;
          for (int r = 0; r < 3; ++r) {
            const double v = direction(r, 0) * dx + direction(r, 1) * dy + direction(r, 2) * dz;
            dist2 += v * v;
          }
          const double dist = std::sqrt(dist2);
          if (!(dist > 0.0) || !std::isfinite(dist)) {
            Rcpp::stop("`vox2ras` is degenerate: two distinct voxels map to the same location.");
          }
          Neighbor nb;
          nb.dx = dx; nb.dy = dy; nb.dz = dz;
          nb.delta = static_cast<std::ptrdiff_t>(dx)
            + static_cast<std::ptrdiff_t>(dy) * static_cast<std::ptrdiff_t>(nx)
            + static_cast<std::ptrdiff_t>(dz) * static_cast<std::ptrdiff_t>(nxy);
          nb.weight = 1.0 / dist;
          nbrs.push_back(nb);
        }
      }
    }
    // The neighborhood sums S_k are bounded by the total neighbor weight, so
    // beta * S_k stays finite (and the log-sum-exp of the MRF term with it)
    // as long as beta times that total is comfortably below the largest
    // double; beyond that the posteriors would degenerate to NaN.
    double total_weight = 0.0;
    for (usize t = 0; t < nbrs.size(); ++t) total_weight += nbrs[t].weight;
    if (!(mrf_beta * total_weight <= 0.5 * std::numeric_limits<double>::max())) {
      Rcpp::stop("`mrf_beta` (%g) is too large for this neighborhood: `mrf_beta` times the total neighbor weight (%g) must stay below %g.",
                 mrf_beta, total_weight, 0.5 * std::numeric_limits<double>::max());
    }
  }

  // ---- EM iterations -------------------------------------------------------
  // The first E-step reads the one-hot initial labels (unlabeled voxels stay
  // all-zero, so they contribute nothing to their neighbors' MRF sums)
  std::vector<double> post_a(N * Ku, 0.0), post_b(N * Ku, 0.0);
  for (usize m = 0; m < N; ++m) {
    if (label[m] >= 0) post_a[m * Ku + static_cast<usize>(label[m])] = 1.0;
  }
  double* post_old = post_a.data();
  double* post_new = post_b.data();

  const usize n_blocks = (N + kBlock - 1) / kBlock;
  const usize n_sums = 4 * Ku + 1;
  std::vector<double> block_sums(n_blocks * n_sums, 0.0);
  std::vector<double> log_norm(Ku), half_inv_var(Ku), log_pi(Ku);
  std::vector<double> trace;
  // `iterations` is only an upper bound (the convergence rule usually stops
  // much earlier), so reserve a small capacity and let the vector grow
  trace.reserve(std::min<usize>(static_cast<usize>(iterations), 64));
  const double two_pi = 6.283185307179586476925286766559;

  if (verbose) {
    Rcpp::Rcout << "[GMM] " << N << " voxels, " << K << " classes, "
                << (use_priors ? "prior" : "k-means") << " initialization (prior weight "
                << w << "), " << nbrs.size() << " MRF neighbors (beta = " << mrf_beta << ")"
                << std::endl;
    Rcpp::Rcout << "[GMM] initial means:";
    for (usize k = 0; k < Ku; ++k) Rcpp::Rcout << " " << mu[k];
    Rcpp::Rcout << std::endl;
  }

  double mmp_prev = -std::numeric_limits<double>::infinity();
  for (int it = 0; it < iterations; ++it) {
    for (usize k = 0; k < Ku; ++k) {
      log_norm[k] = -0.5 * std::log(two_pi * var[k]);
      half_inv_var[k] = 0.5 / var[k];
      log_pi[k] = std::log(std::max(pi[k], kProbEps));
    }

    GmmEStepWorker worker;
    worker.y = y.data();
    worker.grid = grid.data();
    worker.map = map.data();
    worker.nx = nx; worker.ny = ny; worker.nz = nz; worker.nxy = nxy;
    worker.nbrs = &nbrs;
    worker.post_old = post_old;
    worker.prior = prior.empty() ? nullptr : prior.data();
    worker.mu = mu.data();
    worker.log_norm = log_norm.data();
    worker.half_inv_var = half_inv_var.data();
    worker.pi = pi.data();
    worker.log_pi = log_pi.data();
    worker.K = K; worker.N = N; worker.beta = mrf_beta; worker.w = w;
    worker.post_new = post_new;
    worker.block_sums = block_sums.data();
    TinyParallel::parallelFor(0, n_blocks, worker, 1);

    // Combine the block partial sums in a fixed order (thread-independent)
    std::vector<double> s0(Ku, 0.0), s1(Ku, 0.0), s2(Ku, 0.0), sq(Ku, 0.0);
    double max_sum = 0.0;
    for (usize b = 0; b < n_blocks; ++b) {
      const double* bs = block_sums.data() + b * n_sums;
      for (usize k = 0; k < Ku; ++k) {
        s0[k] += bs[4 * k];
        s1[k] += bs[4 * k + 1];
        s2[k] += bs[4 * k + 2];
        sq[k] += bs[4 * k + 3];
      }
      max_sum += bs[4 * Ku];
    }

    // M-step (the sums were taken around the previous means)
    for (usize k = 0; k < Ku; ++k) {
      if (s0[k] > 0.0) {
        const double e1 = s1[k] / s0[k];
        mu[k] += e1;
        var[k] = weighted_variance(s2[k] - e1 * s1[k], s0[k], sq[k], var_floor);
      }
      pi[k] = s0[k] / static_cast<double>(N);
    }

    const double mmp = max_sum / static_cast<double>(N);
    trace.push_back(mmp);
    std::swap(post_old, post_new);

    if (verbose) {
      Rcpp::Rcout << "[GMM] iteration " << (it + 1) << "/" << iterations
                  << ": mean max posterior = " << std::setprecision(6) << mmp
                  << ", means:";
      for (usize k = 0; k < Ku; ++k) Rcpp::Rcout << " " << mu[k];
      Rcpp::Rcout << std::endl;
    }
    Rcpp::checkUserInterrupt();

    // Convergence rule of Atropos: the measure of this iteration is compared
    // with the previous one and the loop stops as soon as the signed change
    // falls below the threshold, so with a zero threshold the first decrease
    // stops it. The first iteration never stops (Atropos compares it against
    // the most negative representable value). The posteriors and labels of
    // the iteration just completed are returned, as in Atropos (no rollback).
    if (it > 0 && (mmp - mmp_prev) < tolerance) {
      if (verbose) {
        Rcpp::Rcout << "[GMM] converged: change " << (mmp - mmp_prev)
                    << " below tolerance " << tolerance << std::endl;
      }
      break;
    }
    mmp_prev = mmp;
  }

  // ---- outputs (post_old holds the latest posteriors) ----------------------
  Rcpp::IntegerVector seg(static_cast<R_xlen_t>(n_full));
  Rcpp::List post_list(K);
  std::vector<double*> post_ptr(Ku);
  for (usize k = 0; k < Ku; ++k) {
    Rcpp::NumericVector p(static_cast<R_xlen_t>(n_full));
    post_list[k] = p;
    post_ptr[k] = REAL(p);
  }
  for (usize m = 0; m < N; ++m) {
    const usize o = grid[m];
    const double* pm = post_old + m * Ku;
    usize best = 0;
    for (usize k = 0; k < Ku; ++k) {
      post_ptr[k][o] = pm[k];
      if (pm[k] > pm[best]) best = k;
    }
    seg[o] = static_cast<int>(best + 1);
  }
  Rcpp::NumericVector means(K), sds(K), props(K);
  for (usize k = 0; k < Ku; ++k) {
    means[k] = mu[k];
    sds[k] = std::sqrt(var[k]);
    props[k] = pi[k];
  }

  return Rcpp::List::create(
    Rcpp::Named("segmentation") = seg,
    Rcpp::Named("posteriors") = post_list,
    Rcpp::Named("means") = means,
    Rcpp::Named("sds") = sds,
    Rcpp::Named("proportions") = props,
    Rcpp::Named("trace") = Rcpp::wrap(trace));
}
