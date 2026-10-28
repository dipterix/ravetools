// mris_outer_surface.cpp
//
// Volume and topology kernels behind mris_outer_surface(), a native port of
// the outer (envelope) surface that the local-gyrification script
// `mris_compute_lgi` writes as `?h.pial-outer-smoothed`:
//
//   mrisFillInteriorParity  fills the lattice points inside a closed surface
//                           (the "voxel center inside" rule of mris_fill) by
//                           scan-line parity;
//   mrisCloseBinary         exact binary morphological closing with a digital
//                           disk (within slices) or ball, computed from
//                           squared Euclidean distance transforms;
//   mrisLargestComponent    largest vertex-connected component of a mesh
//                           (mris_extract_main_component).
//
// Reference: Schaer M, Bach Cuadra M, Tamarit L, Lazeyras F, Eliez S,
// Thiran JP (2008). A surface-based approach to quantify local cortical
// gyrification. IEEE Transactions on Medical Imaging, 27(2), 161-170.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace {

typedef std::int64_t i64;

const double DIST_INF = std::numeric_limits<double>::infinity();

// ---------------------------------------------------------------------------
// One-dimensional squared Euclidean distance transform (lower envelope of
// parabolas; Felzenszwalb and Huttenlocher, 2012):
//
//   d[q] = min_p (q - p)^2 + f[p]
//
// Entries with f = DIST_INF are not sites; a line without sites stays at
// DIST_INF. `v` and `z` are work buffers of length n and n + 1.
// ---------------------------------------------------------------------------
void edt1d(const double *f, double *d, int n, int *v, double *z)
{
    int k = -1;
    for (int q = 0; q < n; q++) {
        if (f[q] == DIST_INF) continue;
        if (k < 0) {
            k = 0; v[0] = q; z[0] = -DIST_INF; z[1] = DIST_INF;
            continue;
        }
        double s;
        for (;;) {
            const int p = v[k];
            s = ((f[q] + (double)q * q) - (f[p] + (double)p * p)) / (2.0 * (q - p));
            // z[0] is -inf, so the search always stops at k = 0 at the latest
            if (s > z[k]) break;
            k--;
        }
        k++;
        v[k] = q; z[k] = s; z[k + 1] = DIST_INF;
    }
    if (k < 0) {
        std::fill(d, d + n, DIST_INF);
        return;
    }
    int j = 0;
    for (int q = 0; q < n; q++) {
        while (z[j + 1] < q) j++;
        const double dq = (double)(q - v[j]);
        d[q] = dq * dq + f[v[j]];
    }
}

i64 axisStride(const int dims[3], int axis)
{
    if (axis == 0) return 1;
    if (axis == 1) return (i64)dims[0];
    return (i64)dims[0] * dims[1];
}

// Apply edt1d to every line of the column-major array `a` along `axis`.
void edtAlongAxis(std::vector<double> &a, const int dims[3], int axis)
{
    const int n = dims[axis];
    const int axis_a = (axis == 0) ? 1 : 0;
    const int axis_b = (axis == 2) ? 1 : 2;
    const i64 stride = axisStride(dims, axis);
    const i64 stride_a = axisStride(dims, axis_a);
    const i64 stride_b = axisStride(dims, axis_b);

    std::vector<double> f(n), d(n), z(n + 1);
    std::vector<int> v(n);

    for (int ib = 0; ib < dims[axis_b]; ib++) {
        for (int ia = 0; ia < dims[axis_a]; ia++) {
            const i64 base = ia * stride_a + ib * stride_b;
            bool has_site = false;
            for (int q = 0; q < n; q++) {
                f[q] = a[base + q * stride];
                if (f[q] != DIST_INF) has_site = true;
            }
            // a line without sites is already all DIST_INF
            if (!has_site) continue;
            edt1d(f.data(), d.data(), n, v.data(), z.data());
            for (int q = 0; q < n; q++) a[base + q * stride] = d[q];
        }
    }
}

// Squared distance from every lattice point to the nearest site, measured
// along every axis except `skip_axis` (-1: all three axes). Leaving one axis
// out makes the transform two-dimensional and independent in each slice
// perpendicular to that axis. Points outside the array are never sites.
std::vector<double> squaredDistanceToSites(const std::vector<unsigned char> &site,
                                           const int dims[3], int skip_axis)
{
    const i64 n = (i64)dims[0] * dims[1] * dims[2];
    std::vector<double> a(n);
    for (i64 i = 0; i < n; i++) a[i] = site[i] ? 0.0 : DIST_INF;
    for (int axis = 0; axis < 3; axis++) {
        if (axis == skip_axis) continue;
        edtAlongAxis(a, dims, axis);
    }
    return a;
}

void readDims(const Rcpp::IntegerVector &dims_, int dims[3], const char *caller)
{
    if (dims_.size() != 3) {
        Rcpp::stop("`%s`: `dims` must have length 3", caller);
    }
    for (int i = 0; i < 3; i++) {
        if (dims_[i] == NA_INTEGER || dims_[i] < 1) {
            Rcpp::stop("`%s`: `dims` must be positive integers", caller);
        }
        dims[i] = dims_[i];
    }
}

} // namespace

// ---------------------------------------------------------------------------
// Fill the lattice points that lie inside a closed triangle surface.
//
// `vb_` holds the vertex positions in continuous lattice-index coordinates
// (lattice point (i, j, k) sits at (i, j, k), zero-based); `it_` is the
// 3 x nf one-based face matrix. Each face is projected onto the i-j plane,
// the columns (i, j) it covers are found with barycentric coordinates, and the
// k-depth of every crossing is recorded; lattice points between consecutive
// pairs of sorted crossings are inside. The columns are tested at
// (i + 1.2345e-4, j + 2.3456e-4): a fixed tiny offset keeps them off the
// projected edges and vertices, so a column crosses a closed surface an even
// number of times. Columns with an odd count (open or self-intersecting
// input) are left empty and counted in the "odd_columns" attribute.
//
// Returns a logical array with dimension `dims`.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::LogicalVector mrisFillInteriorParity(const Rcpp::NumericMatrix &vb_,
                                           const Rcpp::IntegerMatrix &it_,
                                           const Rcpp::IntegerVector &dims)
{
    int d[3];
    readDims(dims, d, "mrisFillInteriorParity");
    if (vb_.nrow() != 3 || it_.nrow() != 3) {
        Rcpp::stop("`mrisFillInteriorParity`: `vb_` and `it_` must have 3 rows");
    }
    const int nv = vb_.ncol(), nf = it_.ncol();
    const int ni = d[0], nj = d[1], nk = d[2];
    const double eps_i = 1.2345e-4, eps_j = 2.3456e-4;

    // (column index i + j * ni, crossing depth along k)
    std::vector<std::pair<i64, double>> hits;
    hits.reserve((size_t)nf);

    for (int fi = 0; fi < nf; fi++) {
        if (fi % 65536 == 0) Rcpp::checkUserInterrupt();
        int a = it_(0, fi) - 1, b = it_(1, fi) - 1, c = it_(2, fi) - 1;
        if (a < 0 || b < 0 || c < 0 || a >= nv || b >= nv || c >= nv) {
            Rcpp::stop("`mrisFillInteriorParity`: face index out of range");
        }
        const double x1 = vb_(0, a) - eps_i, y1 = vb_(1, a) - eps_j, z1 = vb_(2, a);
        const double x2 = vb_(0, b) - eps_i, y2 = vb_(1, b) - eps_j, z2 = vb_(2, b);
        const double x3 = vb_(0, c) - eps_i, y3 = vb_(1, c) - eps_j, z3 = vb_(2, c);

        const double imin = std::max(0.0, std::ceil(std::min(x1, std::min(x2, x3))));
        const double imax = std::min(ni - 1.0, std::floor(std::max(x1, std::max(x2, x3))));
        const double jmin = std::max(0.0, std::ceil(std::min(y1, std::min(y2, y3))));
        const double jmax = std::min(nj - 1.0, std::floor(std::max(y1, std::max(y2, y3))));
        if (imin > imax || jmin > jmax) continue;

        const double det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
        if (det == 0.0) continue;

        for (int j = (int)jmin; j <= (int)jmax; j++) {
            for (int i = (int)imin; i <= (int)imax; i++) {
                const double l1 = ((y2 - y3) * (i - x3) + (x3 - x2) * (j - y3)) / det;
                const double l2 = ((y3 - y1) * (i - x3) + (x1 - x3) * (j - y3)) / det;
                const double l3 = 1.0 - l1 - l2;
                if (l1 >= 0.0 && l2 >= 0.0 && l3 >= 0.0) {
                    hits.emplace_back((i64)i + (i64)j * ni, l1 * z1 + l2 * z2 + l3 * z3);
                }
            }
        }
    }

    std::sort(hits.begin(), hits.end());

    const i64 slice = (i64)ni * nj;
    Rcpp::LogicalVector out((R_xlen_t)(slice * nk));
    int odd_columns = 0;

    for (size_t s = 0; s < hits.size(); ) {
        size_t e = s;
        while (e < hits.size() && hits[e].first == hits[s].first) e++;
        if ((e - s) % 2 == 1) {
            odd_columns++;
            s = e;
            continue;
        }
        const i64 column = hits[s].first;
        for (size_t p = s; p < e; p += 2) {
            const double k0 = std::max(0.0, std::ceil(hits[p].second));
            const double k1 = std::min(nk - 1.0, std::floor(hits[p + 1].second));
            for (int k = (int)k0; k <= (int)k1; k++) {
                out[(R_xlen_t)(column + (i64)k * slice)] = TRUE;
            }
        }
        s = e;
    }

    out.attr("dim") = Rcpp::IntegerVector::create(ni, nj, nk);
    out.attr("odd_columns") = odd_columns;
    return out;
}

// ---------------------------------------------------------------------------
// Binary morphological closing (dilation followed by erosion) with the
// digital disk or ball {offset : |offset|^2 <= radius^2}.
//
// The dilation keeps the points within `radius` of the mask; the erosion then
// keeps the points farther than `radius` from every point outside the dilated
// mask. Points outside the array count as outside for the dilation and as
// inside for the erosion (the padding MATLAB's imclose uses).
//
// `slice_axis` (zero-based) names the axis the element does not extend
// along, giving a two-dimensional closing applied independently in every
// slice perpendicular to it; -1 gives the three-dimensional ball.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::LogicalVector mrisCloseBinary(const Rcpp::LogicalVector &mask,
                                    const Rcpp::IntegerVector &dims,
                                    double radius, int slice_axis)
{
    int d[3];
    readDims(dims, d, "mrisCloseBinary");
    const i64 n = (i64)d[0] * d[1] * d[2];
    if ((i64)mask.size() != n) {
        Rcpp::stop("`mrisCloseBinary`: `mask` length does not match `dims`");
    }
    if (!(radius >= 0.0) || !std::isfinite(radius)) {
        Rcpp::stop("`mrisCloseBinary`: `radius` must be a non-negative number");
    }
    if (slice_axis < -1 || slice_axis > 2) {
        Rcpp::stop("`mrisCloseBinary`: `slice_axis` must be -1, 0, 1, or 2");
    }
    const double r2 = radius * radius;

    std::vector<unsigned char> site(n);
    for (i64 i = 0; i < n; i++) site[i] = (mask[(R_xlen_t)i] == TRUE) ? 1 : 0;

    // dilation; afterwards `site` marks the points outside the dilated mask
    std::vector<double> dist = squaredDistanceToSites(site, d, slice_axis);
    for (i64 i = 0; i < n; i++) site[i] = (dist[i] <= r2) ? 0 : 1;
    Rcpp::checkUserInterrupt();

    // erosion of the dilated mask
    dist = squaredDistanceToSites(site, d, slice_axis);
    Rcpp::LogicalVector out((R_xlen_t)n);
    for (i64 i = 0; i < n; i++) out[(R_xlen_t)i] = dist[i] > r2;

    out.attr("dim") = Rcpp::IntegerVector::create(d[0], d[1], d[2]);
    return out;
}

// ---------------------------------------------------------------------------
// Largest vertex-connected component of a triangle mesh, as
// mris_extract_main_component keeps it.
//
// `it_` is the 3 x nf one-based face matrix and `nv` the vertex count.
// Returns a list with `keep` (logical, per vertex) and `n_components` (the
// number of connected components, isolated vertices included). Ties go to
// the component that holds the lowest vertex index.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mrisLargestComponent(const Rcpp::IntegerMatrix &it_, int nv)
{
    if (it_.nrow() != 3) {
        Rcpp::stop("`mrisLargestComponent`: `it_` must have 3 rows");
    }
    if (nv < 1) {
        Rcpp::stop("`mrisLargestComponent`: `nv` must be positive");
    }
    const int nf = it_.ncol();

    // union-find whose roots are the lowest vertex index of each component
    std::vector<int> root(nv);
    std::iota(root.begin(), root.end(), 0);
    auto find = [&root](int i) {
        while (root[i] != i) {
            root[i] = root[root[i]];
            i = root[i];
        }
        return i;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if (a < b) root[b] = a; else if (b < a) root[a] = b;
    };

    for (int fi = 0; fi < nf; fi++) {
        const int a = it_(0, fi) - 1, b = it_(1, fi) - 1, c = it_(2, fi) - 1;
        if (a < 0 || b < 0 || c < 0 || a >= nv || b >= nv || c >= nv) {
            Rcpp::stop("`mrisLargestComponent`: face index out of range");
        }
        unite(a, b);
        unite(b, c);
    }

    std::vector<int> size(nv, 0);
    for (int v = 0; v < nv; v++) size[find(v)]++;

    int n_components = 0, best = -1;
    for (int v = 0; v < nv; v++) {
        if (root[v] != v) continue;
        n_components++;
        if (best < 0 || size[v] > size[best]) best = v;
    }

    Rcpp::LogicalVector keep(nv);
    for (int v = 0; v < nv; v++) keep[v] = (find(v) == best);

    return Rcpp::List::create(
        Rcpp::Named("keep") = keep,
        Rcpp::Named("n_components") = n_components
    );
}
