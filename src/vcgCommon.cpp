#include <Rcpp.h>
#include "vcgCommon.h"
#include <vcg/complex/algorithms/hole.h>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <vector>

namespace {

// Count boundary (multiplicity 1) and non-manifold (multiplicity > 2) edges
// directly from the face/vertex-index structure, robust regardless of
// whether VCG's FaceFace topology has been (re)built. Shared by routines
// that need to detect or report mesh topology defects (e.g. vcgFixDefects,
// vcgCountEdgeDefects).
void countEdgeDefects(ravetools::MyMesh &m, int &boundary, int &nonmanifold)
{
  std::unordered_map<int64_t, int> edge_count;
  edge_count.reserve((size_t)m.fn * 3);
  int64_t nv = (int64_t)m.vert.size();

  auto add_edge = [&](int a, int b) {
    if (a > b) std::swap(a, b);
    edge_count[(int64_t)a * nv + (int64_t)b]++;
  };

  for (ravetools::MyMesh::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
    if (fi->IsD()) continue;
    int i0 = (int)vcg::tri::Index(m, fi->V(0));
    int i1 = (int)vcg::tri::Index(m, fi->V(1));
    int i2 = (int)vcg::tri::Index(m, fi->V(2));
    add_edge(i0, i1);
    add_edge(i1, i2);
    add_edge(i2, i0);
  }

  boundary = 0;
  nonmanifold = 0;
  for (const auto &kv : edge_count) {
    if (kv.second == 1) boundary++;
    else if (kv.second > 2) nonmanifold++;
  }
}

// Average length of all (non-deleted) face edges of the mesh, counting each
// edge once per incident face (i.e. shared edges are counted twice). Useful
// as a scale-aware reference length, e.g. to derive a vertex-merge tolerance.
double averageEdgeLength(ravetools::MyMesh &m)
{
  double sum = 0.0;
  int cnt = 0;
  for (ravetools::MyMesh::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
    if (fi->IsD()) continue;
    for (int j = 0; j < 3; j++) {
      sum += vcg::Distance(fi->V(j)->P(), fi->V((j + 1) % 3)->P());
      cnt++;
    }
  }
  return (cnt > 0) ? (sum / cnt) : 0.0;
}

double maxEdgeLength(ravetools::MyMesh &m)
{
  double maxLen = 0.0;
  for (ravetools::MyMesh::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
    if (fi->IsD()) continue;
    for (int j = 0; j < 3; j++) {
      double len = vcg::Distance(fi->V(j)->P(), fi->V((j + 1) % 3)->P());
      if (len > maxLen) maxLen = len;
    }
  }
  return maxLen;
}

// ---- implicit smoothing (vcgSmoothImplicit) --------------------------------
//
// Smoothing solves, per coordinate axis, (M + lambda * L^k) X = M V with
// k = 2^(degree - 1): M is the per-vertex sum of doubled incident face areas
// scaled by its maximum, and L the Laplacian accumulated per face edge. This
// is the system vcglib's ImplicitSmoother assembles, but vcglib stacks the
// three axes into one 3N x 3N system and factorizes it with 32-bit indices:
// on a 3.2 M vertex whole-brain isosurface (degree 2) the factor needs
// 2.85e9 nonzeros, the index sum wraps negative, the factor is never
// allocated and the factorization writes through a null pointer. Here the
// axes share one N x N matrix, in double precision with 64-bit indices, and
// every allocation proportional to the mesh is estimated before it is made.

typedef Eigen::Index SmoothIndex;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, SmoothIndex> SmoothMatrix;

// Symmetric sparsity pattern in compressed columns, diagonal included and
// row indices sorted
struct SmoothPattern {
  std::vector<int64_t> ptr;
  std::vector<int> idx;
};

double smoothMatrixBytes(double nnz, double n)
{
  return nnz * (sizeof(double) + sizeof(SmoothIndex)) + (n + 1.0) * sizeof(SmoothIndex);
}

double smoothPatternBytes(double nnz, double n)
{
  return nnz * sizeof(int) + (n + 1.0) * sizeof(int64_t);
}

// Peak memory of the solve: each stage records the bytes it holds at once
// (live arrays, not what the allocator may keep after they are freed), and
// `check` stops with an explanation when the largest exceeds the limit
struct SmoothMemoryGuard {
  double limit;
  int nVertices;
  int degree;
  double base;   // the input arrays, held throughout
  double peak;

  void need(double bytes) {
    if (base + bytes > peak) peak = base + bytes;
  }

  // `partial`: stages still to come are not counted yet
  void check(bool partial) {
    if (peak <= limit) return;
    const double gib = 1024.0 * 1024.0 * 1024.0;
    Rcpp::stop(
      "vcg_smooth_implicit: smoothing this mesh (%d vertices, degree %d) "
      "needs %s %.3g GiB of memory, more than `max_memory` (%.3g GiB). "
      "Reduce the mesh first (for example with `vcg_decimate()`), use a "
      "lower `degree`, smooth explicitly with `vcg_smooth_explicit()` or "
      "`mris_smooth()`, or raise `max_memory` if this machine has the memory.",
      nVertices, degree, partial ? "at least" : "about", peak / gib, limit / gib);
  }
};

// Number of nonzeros of A * A for a symmetric pattern A. When `out` is not
// null the (sorted) pattern of the product is stored there as well.
int64_t squarePattern(const SmoothPattern &a, int n, SmoothPattern *out)
{
  std::vector<int> mark(n, -1);
  int64_t total = 0;
  if (out) out->ptr.assign((size_t)n + 1, 0);
  for (int c = 0; c < n; c++) {
    for (int64_t k = a.ptr[c]; k < a.ptr[c + 1]; k++) {
      const int j = a.idx[k];
      for (int64_t l = a.ptr[j]; l < a.ptr[j + 1]; l++) {
        const int i = a.idx[l];
        if (mark[i] != c) { mark[i] = c; total++; }
      }
    }
    if (out) out->ptr[c + 1] = total;
  }
  if (!out) return total;

  out->idx.resize((size_t)total);
  std::fill(mark.begin(), mark.end(), -1);
  for (int c = 0; c < n; c++) {
    int64_t pos = out->ptr[c];
    for (int64_t k = a.ptr[c]; k < a.ptr[c + 1]; k++) {
      const int j = a.idx[k];
      for (int64_t l = a.ptr[j]; l < a.ptr[j + 1]; l++) {
        const int i = a.idx[l];
        if (mark[i] != c) { mark[i] = c; out->idx[pos++] = i; }
      }
    }
    std::sort(out->idx.begin() + out->ptr[c], out->idx.begin() + pos);
  }
  return total;
}

// A * A for a symmetric matrix with sorted columns (Gustavson's algorithm
// with a dense accumulator), producing sorted columns. Eigen's own sparse
// product keeps three copies of the result alive while it sorts them.
SmoothMatrix squareSymmetric(const SmoothMatrix &a)
{
  const SmoothIndex n = a.cols();
  const SmoothIndex *ap = a.outerIndexPtr();
  const SmoothIndex *ai = a.innerIndexPtr();
  const double *ax = a.valuePtr();

  std::vector<SmoothIndex> mark(n, -1);
  std::vector<SmoothIndex> colPtr(n + 1, 0);
  for (SmoothIndex c = 0; c < n; c++) {
    SmoothIndex cnt = 0;
    for (SmoothIndex k = ap[c]; k < ap[c + 1]; k++) {
      const SmoothIndex j = ai[k];
      for (SmoothIndex l = ap[j]; l < ap[j + 1]; l++) {
        const SmoothIndex i = ai[l];
        if (mark[i] != c) { mark[i] = c; cnt++; }
      }
    }
    colPtr[c + 1] = colPtr[c] + cnt;
  }

  SmoothMatrix r(n, n);
  r.resizeNonZeros(colPtr[n]);
  SmoothIndex *rp = r.outerIndexPtr();
  SmoothIndex *ri = r.innerIndexPtr();
  double *rx = r.valuePtr();
  std::copy(colPtr.begin(), colPtr.end(), rp);
  std::vector<SmoothIndex>().swap(colPtr);

  std::fill(mark.begin(), mark.end(), -1);
  std::vector<double> acc(n, 0.0);
  for (SmoothIndex c = 0; c < n; c++) {
    SmoothIndex pos = rp[c];
    for (SmoothIndex k = ap[c]; k < ap[c + 1]; k++) {
      const SmoothIndex j = ai[k];
      const double ajc = ax[k];
      for (SmoothIndex l = ap[j]; l < ap[j + 1]; l++) {
        const SmoothIndex i = ai[l];
        if (mark[i] != c) { mark[i] = c; ri[pos++] = i; acc[i] = 0.0; }
        acc[i] += ax[l] * ajc;
      }
    }
    std::sort(ri + rp[c], ri + rp[c + 1]);
    for (SmoothIndex k = rp[c]; k < rp[c + 1]; k++) rx[k] = acc[ri[k]];
  }
  return r;
}

// Reference to the stored entry (row, col); the entry must exist
double &smoothEntry(SmoothMatrix &m, SmoothIndex row, SmoothIndex col)
{
  SmoothIndex *begin = m.innerIndexPtr() + m.outerIndexPtr()[col];
  SmoothIndex *end = m.innerIndexPtr() + m.outerIndexPtr()[col + 1];
  SmoothIndex *hit = std::lower_bound(begin, end, row);
  return m.valuePtr()[hit - m.innerIndexPtr()];
}

// Nonzeros of the strictly lower LDLT factor of `ap` (upper triangle stored,
// already permuted), counted in 64 bits with the elimination-tree walk of
// Eigen's SimplicialCholeskyBase::analyzePattern_preordered, which sums the
// same counts in the matrix's index type
double ldltFactorNonZeros(const SmoothMatrix &ap)
{
  const SmoothIndex n = ap.cols();
  std::vector<SmoothIndex> parent(n), tags(n);
  double total = 0.0;
  for (SmoothIndex k = 0; k < n; ++k) {
    parent[k] = -1;
    tags[k] = k;
    for (SmoothMatrix::InnerIterator it(ap, k); it; ++it) {
      SmoothIndex i = it.index();
      if (i < k) {
        for (; tags[i] != k; i = parent[i]) {
          if (parent[i] == -1) parent[i] = k;
          total += 1.0;
          tags[i] = k;
        }
      }
    }
  }
  return total;
}

int findRoot(std::vector<int> &root, int i)
{
  while (root[i] != i) {
    root[i] = root[root[i]];
    i = root[i];
  }
  return i;
}

} // namespace

// [[Rcpp::export]]
SEXP vcgIsoSurface(SEXP array_, double thresh) {
  try {
    Rcpp::IntegerVector arrayDims( Rf_getAttrib(array_, R_DimSymbol) );
    std::vector<float> vecArray = Rcpp::as<std::vector<float> >(array_);

    ravetools::MyMesh m;
    ravetools::VertexIterator vi;
    ravetools::FaceIterator fi;
    int i,j,k;

    // typedef MySimpleVolume<ravetools::MySimpleVoxel> MyVolume;
    ravetools::MyVolume	volume;
    // typedef vcg::tri::TrivialWalker<ravetools::MyMesh, ravetools::MyVolume>	MyWalker;
    // typedef vcg::tri::MarchingCubes<ravetools::MyMesh, ravetools::MyWalker>	MyMarchingCubes;
    ravetools::MyWalker walker;
    volume.Init( vcg::Point3i(arrayDims[0], arrayDims[1], arrayDims[2]) );
    for( i = 0; i < arrayDims[0]; i++ ) {
      for( j = 0; j < arrayDims[1]; j++ ) {
        for( k = 0; k < arrayDims[2]; k++ ) {
          int tmpval = vecArray[ i + j * arrayDims[0] + k * ( arrayDims[0] * arrayDims[1] ) ];
          /*if (tmpval >= lower && tmpval <= upper)
           volume.Val(i,j,k)=tmpval;
           else*/
          volume.Val(i,j,k)=tmpval;
        }
      }
    }

    Rcpp::checkUserInterrupt();
    //write back
    /*volume.Init(Point3i(64,64,64));
     for(int i=0;i<64;i++)
     for(int j=0;j<64;j++)
     for(int k=0;k<64;k++)
     volume.Val(i,j,k)=(j-32)*(j-32)+(k-32)*(k-32)  + i*10*(float)math::Perlin::Noise(i*.2,j*.2,k*.2);*/
    ravetools::MyMarchingCubes	mc(m, walker);
    walker.BuildMesh<ravetools::MyMarchingCubes>(m, volume, mc, thresh);
    vcg::tri::Allocator< ravetools::MyMesh >::CompactVertexVector(m);
    vcg::tri::Allocator< ravetools::MyMesh >::CompactFaceVector(m);
    vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal< ravetools::MyMesh >::NormalizePerVertex(m);
    vcg::SimpleTempData< ravetools::MyMesh::VertContainer, int > indiceout(m.vert);
    Rcpp::NumericMatrix vbout(3,m.vn), normals(3,m.vn);
    Rcpp::IntegerMatrix itout(3,m.fn);

    Rcpp::checkUserInterrupt();

    vi=m.vert.begin();
    for (i=0;  i < m.vn; i++) {
      indiceout[vi] = i;
      vbout(0,i) = (*vi).P()[0];
      vbout(1,i) = (*vi).P()[1];
      vbout(2,i) = (*vi).P()[2];
      normals(0,i) = (*vi).N()[0];
      normals(1,i) = (*vi).N()[1];
      normals(2,i) = (*vi).N()[2];
      ++vi;
    }
    ravetools::FacePointer fp;

    fi=m.face.begin();
    j = 0;
    for (i=0; i < m.fn; i++) {
      fp=&(*fi);
      itout(0,i) = indiceout[fp->cV(0)]+1;
      itout(1,i) = indiceout[fp->cV(1)]+1;
      itout(2,i) = indiceout[fp->cV(2)]+1;
      ++fi;
    }
    //delete &walker;
    //delete &volume;
    return Rcpp::List::create(Rcpp::Named("vb") = vbout,
                              Rcpp::Named("it") = itout,
                              Rcpp::Named("normals") = normals);
  } catch (std::exception& e) {
    Rcpp::stop( e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgSmoothImplicit(
    SEXP vb_, SEXP it_, double lambda, bool useMassMatrix, bool fixBorder,
    bool useCotWeight, int degree, double lapWeight, bool SmoothQ,
    double maxMemory)
{
  try {
    if (SmoothQ) {
      Rcpp::stop("vcg_smooth_implicit: smoothing per-vertex quality is not supported");
    }
    if (!std::isfinite(lambda) || lambda < 0.0) {
      Rcpp::stop("vcg_smooth_implicit: `lambda` must be a non-negative number");
    }
    if (degree < 1) degree = 1;

    Rcpp::NumericMatrix vb(vb_);
    Rcpp::IntegerMatrix it(it_);
    const int n = vb.ncol();
    const int nf = it.ncol();
    if (vb.nrow() < 3 || it.nrow() != 3) {
      Rcpp::stop("vcg_smooth_implicit: `vb` must have 3 rows and `it` exactly 3 rows");
    }
    for (R_xlen_t k = 0; k < it.size(); k++) {
      if (it[k] < 0 || it[k] >= n) {
        Rcpp::stop("vcg_smooth_implicit: face index out of range");
      }
    }

    const double dn = n, df = nf, vec = n * (double)sizeof(double);
    SmoothMemoryGuard guard = { maxMemory, n, degree,
                                3.0 * vec + 3.0 * df * sizeof(int), 0.0 };

    // ---- vertex adjacency: the pattern of L, border vertices --------------
    // every face edge (a, b) lists b in column a and a in column b; those
    // slots and the merged pattern (at most n + 6 nf entries) coexist
    guard.need(6.0 * df * sizeof(int) + 2.0 * (dn + 1.0) * sizeof(int64_t) + dn
               + smoothPatternBytes(dn + 6.0 * df, dn));
    guard.check(true);
    SmoothPattern pattern;
    std::vector<char> border(n, 0);
    {
      std::vector<int64_t> slot((size_t)n + 1, 0);
      for (int f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
          const int a = it(e, f), b = it((e + 1) % 3, f);
          if (a == b) continue;
          slot[a + 1]++;
          slot[b + 1]++;
        }
      }
      for (int i = 0; i < n; i++) slot[i + 1] += slot[i];
      std::vector<int> nbr((size_t)slot[n]);
      std::vector<int64_t> cursor(slot.begin(), slot.end() - 1);
      for (int f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
          const int a = it(e, f), b = it((e + 1) % 3, f);
          if (a == b) continue;
          nbr[cursor[a]++] = b;
          nbr[cursor[b]++] = a;
        }
      }
      std::vector<int64_t>().swap(cursor);

      // sort each column, merge repeated neighbors (an edge used by exactly
      // one face lies on the border) and insert the diagonal
      int64_t nnz = 0;
      for (int c = 0; c < n; c++) {
        std::sort(nbr.begin() + slot[c], nbr.begin() + slot[c + 1]);
        int64_t k = slot[c];
        while (k < slot[c + 1]) {
          int64_t run = k;
          while (run < slot[c + 1] && nbr[run] == nbr[k]) run++;
          if (run - k == 1) border[c] = 1;
          nnz++;
          k = run;
        }
        nnz++;
      }
      pattern.ptr.assign((size_t)n + 1, 0);
      pattern.idx.resize((size_t)nnz);
      int64_t pos = 0;
      for (int c = 0; c < n; c++) {
        bool diag = false;
        int64_t k = slot[c];
        while (k < slot[c + 1]) {
          const int j = nbr[k];
          if (!diag && j > c) { pattern.idx[pos++] = c; diag = true; }
          pattern.idx[pos++] = j;
          while (k < slot[c + 1] && nbr[k] == j) k++;
        }
        if (!diag) pattern.idx[pos++] = c;
        pattern.ptr[c + 1] = pos;
      }
    }
    Rcpp::checkUserInterrupt();

    // ---- vertices held at their input positions ---------------------------
    // a vertex on no face has nothing to smooth against
    std::vector<char> pinned(n, 0);
    for (int c = 0; c < n; c++) {
      const bool isolated = pattern.ptr[c + 1] - pattern.ptr[c] == 1;
      pinned[c] = isolated || (fixBorder && border[c]);
    }
    std::vector<char>().swap(border);

    if (!useMassMatrix) {
      // without the data term only the pinned vertices hold the mesh, so every
      // connected component needs at least one of them
      std::vector<int> root(n);
      for (int i = 0; i < n; i++) root[i] = i;
      for (int c = 0; c < n; c++) {
        for (int64_t k = pattern.ptr[c]; k < pattern.ptr[c + 1]; k++) {
          const int a = findRoot(root, c), b = findRoot(root, pattern.idx[k]);
          if (a != b) root[a] = b;
        }
      }
      std::vector<char> held(n, 0), seen(n, 0);
      for (int i = 0; i < n; i++) {
        if (pinned[i]) held[findRoot(root, i)] = 1;
      }
      int loose = 0, parts = 0;
      for (int i = 0; i < n; i++) {
        const int r = findRoot(root, i);
        if (seen[r]) continue;
        seen[r] = 1;
        parts++;
        if (!held[r]) loose++;
      }
      if (loose > 0) {
        Rcpp::stop(
          "vcg_smooth_implicit: with `use_mass_matrix = FALSE` nothing keeps "
          "vertices near their input positions except the fixed border "
          "vertices, but %d of %d connected parts of the mesh have none "
          "(closed parts have no border, and `fix_border = FALSE` fixes "
          "nothing). Use `fix_border = TRUE` on an open mesh, or "
          "`use_mass_matrix = TRUE`.", loose, parts);
      }
    }

    // ---- memory needed, from the exact sizes of the powers of L -----------
    const int squarings = degree - 1;
    const bool direct = degree >= 3;
    std::vector<double> nnzPow(1, (double)pattern.idx.size());
    {
      SmoothPattern cur;
      const SmoothPattern *src = &pattern;
      for (int s = 1; s <= squarings; s++) {
        const double cnt = (double)squarePattern(*src, n, NULL);
        nnzPow.push_back(cnt);
        if (s < squarings) {
          // the next count needs this power's pattern
          guard.need(smoothPatternBytes(nnzPow[0], dn)
                     + (s > 1 ? smoothPatternBytes(nnzPow[s - 1], dn) : 0.0)
                     + smoothPatternBytes(cnt, dn));
          guard.check(true);
          SmoothPattern next;
          squarePattern(*src, n, &next);
          cur.ptr.swap(next.ptr);
          cur.idx.swap(next.idx);
          src = &cur;
        }
        Rcpp::checkUserInterrupt();
      }
    }
    const double nnzK = nnzPow.back();

    // L next to its pattern, then each squaring with its operand
    guard.need(smoothPatternBytes(nnzPow[0], dn) + smoothMatrixBytes(nnzPow[0], dn) + vec);
    for (int s = 1; s <= squarings; s++) {
      guard.need(smoothMatrixBytes(nnzPow[s - 1], dn) + smoothMatrixBytes(nnzPow[s], dn)
                 + 3.0 * vec + vec);
    }
    // the system with its right-hand sides, initial guesses and solutions
    const double system = smoothMatrixBytes(nnzK, dn) + 10.0 * vec + dn;
    if (direct) {
      const double nnzUpper = (nnzK - dn) / 2.0 + dn;
      // AMD copies the symmetric pattern and grows it by a fifth
      guard.need(system + 2.2 * nnzK * (sizeof(double) + sizeof(SmoothIndex))
                 + 10.0 * (dn + 1.0) * sizeof(SmoothIndex));
      guard.need(system + smoothMatrixBytes(nnzUpper, dn) + 2.0 * dn * sizeof(SmoothIndex));
    } else {
      guard.need(system + 5.0 * vec);
    }
    // normals are computed on a vcg mesh rebuilt from the result
    guard.need(dn * sizeof(ravetools::MyVertex) + df * sizeof(ravetools::MyFace)
               + dn * (sizeof(void*) + sizeof(unsigned int)) + df * sizeof(unsigned int)
               + 6.0 * vec + 3.0 * df * sizeof(int));
    // the factor of the direct path is counted once its ordering is known
    guard.check(direct);

    // ---- assemble L and its powers ----------------------------------------
    SmoothMatrix S(n, n);
    S.resizeNonZeros(pattern.idx.size());
    std::copy(pattern.ptr.begin(), pattern.ptr.end(), S.outerIndexPtr());
    std::copy(pattern.idx.begin(), pattern.idx.end(), S.innerIndexPtr());
    std::fill(S.valuePtr(), S.valuePtr() + S.nonZeros(), 0.0);
    std::vector<int64_t>().swap(pattern.ptr);
    std::vector<int>().swap(pattern.idx);

    std::vector<double> mass(n, 0.0);
    for (int f = 0; f < nf; f++) {
      const int idx[3] = { it(0, f), it(1, f), it(2, f) };
      double p[3][3];
      for (int v = 0; v < 3; v++) {
        for (int d = 0; d < 3; d++) p[v][d] = vb(d, idx[v]);
      }
      const double u[3] = { p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2] };
      const double w[3] = { p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2] };
      const double cr[3] = { u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                             u[0] * w[1] - u[1] * w[0] };
      const double doubleArea = std::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
      for (int e = 0; e < 3; e++) {
        mass[idx[e]] += doubleArea;
        const int a = idx[e], b = idx[(e + 1) % 3], o = (e + 2) % 3;
        if (a == b) continue;
        double weight = lapWeight;
        if (useCotWeight) {
          // half the cotangent of the angle opposite the edge, as
          // vcg::tri::Harmonic::CotangentWeight returns without
          // face-face adjacency
          double ca[3], cb[3];
          for (int d = 0; d < 3; d++) {
            ca[d] = p[e][d] - p[o][d];
            cb[d] = p[(e + 1) % 3][d] - p[o][d];
          }
          const double dot = ca[0] * cb[0] + ca[1] * cb[1] + ca[2] * cb[2];
          const double cx = ca[1] * cb[2] - ca[2] * cb[1];
          const double cy = ca[2] * cb[0] - ca[0] * cb[2];
          const double cz = ca[0] * cb[1] - ca[1] * cb[0];
          weight = dot / std::sqrt(cx * cx + cy * cy + cz * cz) / 2.0;
        }
        smoothEntry(S, a, a) += weight;
        smoothEntry(S, b, b) += weight;
        smoothEntry(S, a, b) -= weight;
        smoothEntry(S, b, a) -= weight;
      }
    }
    Rcpp::checkUserInterrupt();

    for (int s = 1; s <= squarings; s++) {
      SmoothMatrix next = squareSymmetric(S);
      S.swap(next);
      Rcpp::checkUserInterrupt();
    }

    // ---- S = M + lambda * L^k, pinned vertices as identity rows -----------
    double maxMass = 0.0;
    for (int i = 0; i < n; i++) maxMass = std::max(maxMass, mass[i]);
    if (useMassMatrix && !(maxMass > 0.0)) {
      Rcpp::stop("vcg_smooth_implicit: every face of the mesh has zero area");
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, 3);
    Eigen::MatrixXd X0(n, 3);
    for (int i = 0; i < n; i++) {
      for (int d = 0; d < 3; d++) X0(i, d) = vb(d, i);
    }
    {
      double *sx = S.valuePtr();
      for (SmoothIndex k = 0; k < S.nonZeros(); k++) sx[k] *= lambda;
    }
    if (useMassMatrix) {
      for (int i = 0; i < n; i++) {
        const double m = mass[i] / maxMass;
        smoothEntry(S, i, i) += m;
        B.row(i) = m * X0.row(i);
      }
    }
    std::vector<double>().swap(mass);

    // move each pinned vertex's couplings to the right-hand side
    for (int p = 0; p < n; p++) {
      if (!pinned[p]) continue;
      for (SmoothMatrix::InnerIterator iter(S, p); iter; ++iter) {
        const SmoothIndex i = iter.index();
        if (i == p) continue;
        if (!pinned[i]) B.row(i) -= iter.value() * X0.row(p);
        iter.valueRef() = 0.0;
        smoothEntry(S, p, i) = 0.0;
      }
      smoothEntry(S, p, p) = 1.0;
      B.row(p) = X0.row(p);
    }
    S.prune(0.0);

    {
      const double *sx = S.valuePtr();
      bool finite = B.allFinite();
      for (SmoothIndex k = 0; finite && k < S.nonZeros(); k++) finite = std::isfinite(sx[k]);
      if (!finite) {
        Rcpp::stop(
          "vcg_smooth_implicit: the smoothing system has non-finite entries; "
          "the mesh probably has degenerate (zero-area) faces%s",
          useCotWeight ? ", whose cotangent weights are infinite" : "");
      }
    }
    Rcpp::checkUserInterrupt();

    // ---- solve ---------------------------------------------------------------
    Eigen::MatrixXd X(n, 3);
    if (!direct) {
      // conjugate gradients: memory linear in the mesh, and degree <= 2 is
      // well conditioned enough to converge in a few hundred iterations
      Eigen::ConjugateGradient<SmoothMatrix, Eigen::Lower | Eigen::Upper> cg;
      cg.setTolerance(1e-10);
      cg.setMaxIterations(10000);
      cg.compute(S);
      for (int d = 0; d < 3; d++) {
        X.col(d) = cg.solveWithGuess(B.col(d), X0.col(d));
        if (cg.info() != Eigen::Success) {
          Rcpp::stop(
            "vcg_smooth_implicit: the iterative solve did not converge in %d "
            "iterations (relative residual %.3g); try a smaller `lambda` or "
            "`degree`", (int)cg.iterations(), (double)cg.error());
        }
        Rcpp::checkUserInterrupt();
      }
      SmoothMatrix().swap(S);
    } else {
      // degree >= 3 squares L at least twice, too stiff for unpreconditioned
      // conjugate gradients: factorize, after counting the factor's size
      Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, SmoothIndex> perm, permInv;
      {
        Eigen::AMDOrdering<SmoothIndex> amd;
        amd(S.selfadjointView<Eigen::Lower>(), permInv);
      }
      perm = permInv.inverse();
      SmoothMatrix upper(n, n);
      upper.selfadjointView<Eigen::Upper>() = S.selfadjointView<Eigen::Lower>().twistedBy(perm);
      SmoothMatrix().swap(S);
      Rcpp::checkUserInterrupt();

      const double nnzFactor = ldltFactorNonZeros(upper);
      guard.need(smoothMatrixBytes((double)upper.nonZeros(), dn)
                 + smoothMatrixBytes(nnzFactor, dn) + 8.0 * vec
                 + 2.0 * dn * sizeof(SmoothIndex) + 9.0 * vec);
      guard.check(false);

      Eigen::SimplicialLDLT<SmoothMatrix, Eigen::Upper, Eigen::NaturalOrdering<SmoothIndex> > ldlt;
      ldlt.compute(upper);
      if (ldlt.info() != Eigen::Success) {
        Rcpp::stop("vcg_smooth_implicit: the smoothing system is singular");
      }
      for (int d = 0; d < 3; d++) {
        Eigen::VectorXd rhs = perm * B.col(d);
        Eigen::VectorXd sol = ldlt.solve(rhs);
        X.col(d) = permInv * sol;
        Rcpp::checkUserInterrupt();
      }
    }
    if (!X.allFinite()) {
      Rcpp::stop("vcg_smooth_implicit: the solution has non-finite coordinates");
    }

    // pinned vertices keep their input coordinates exactly
    Rcpp::NumericMatrix vbout(3, n);
    for (int i = 0; i < n; i++) {
      for (int d = 0; d < 3; d++) vbout(d, i) = pinned[i] ? vb(d, i) : X(i, d);
    }
    X.resize(0, 0);
    B.resize(0, 0);
    X0.resize(0, 0);

    // ---- normals of the smoothed mesh -------------------------------------
    ravetools::MyMesh m;
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vbout, it_);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(m);
    Rcpp::NumericMatrix normals(3, n);
    for (int i = 0; i < n; i++) {
      for (int d = 0; d < 3; d++) normals(d, i) = m.vert[i].N()[d];
    }
    Rcpp::IntegerMatrix itout(3, nf);
    for (R_xlen_t k = 0; k < itout.size(); k++) itout[k] = it[k] + 1;

    return Rcpp::List::create(Rcpp::Named("vb") = vbout,
                              Rcpp::Named("normals") = normals,
                              Rcpp::Named("it") = itout
    );

  } catch (std::exception& e) {
    Rcpp::stop( e.what());
    return Rcpp::wrap(1);
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}


// [[Rcpp::export]]
SEXP vcgSmooth(SEXP vb_, SEXP it_, int iter, int method, float lambda, float mu, float delta_)
{
  try {
    int i;
    ravetools::MyMesh m;
    ravetools::VertexIterator vi;
    ravetools::FaceIterator fi;
    //set up parameters
    ravetools::ScalarType delta = delta_;
    //allocate mesh and fill it
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m,vb_,it_);

    Rcpp::checkUserInterrupt();

    if (method == 0) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      unsigned int cnt = vcg::tri::UpdateSelection<ravetools::MyMesh>::VertexFromFaceStrict(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordTaubin(m, iter, lambda, mu, cnt>0);
    } else if (method == 1) {
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordLaplacian(m, iter);
    } else if (method == 2) {
      vcg::tri::UpdateSelection<ravetools::MyMesh>::FaceAll(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      unsigned int cnt=vcg::tri::UpdateSelection<ravetools::MyMesh>::VertexFromFaceStrict(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordLaplacianHC(m, iter,cnt>0);
    } else if (method == 3) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceClearB(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordScaleDependentLaplacian_Fujiwara(m,iter,delta);
    } else if (method == 4) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceClearB(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordLaplacianAngleWeighted(m,iter,delta);
    }
    else if (method == 5) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceClearB(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordPlanarLaplacian(m, iter, delta);
    }

    Rcpp::checkUserInterrupt();

    vcg::tri::Allocator<ravetools::MyMesh>::CompactVertexVector(m);
    vcg::tri::Allocator<ravetools::MyMesh>::CompactFaceVector(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(m);
    Rcpp::NumericMatrix vb(3, m.vn);
    Rcpp::NumericMatrix normals(3, m.vn);
    Rcpp::IntegerMatrix itout(3, m.fn);
    //write back output
    vcg::SimpleTempData<ravetools::MyMesh::VertContainer,int>indices(m.vert);

    // write back updated mesh
    vi=m.vert.begin();
    for (i=0; i < m.vn; i++) {
      indices[vi] = i;
      if( ! vi->IsD() ) {
        vb(0,i) = (*vi).P()[0];
        vb(1,i) = (*vi).P()[1];
        vb(2,i) = (*vi).P()[2];
        normals(0,i) = (*vi).N()[0];
        normals(1,i) = (*vi).N()[1];
        normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }

    ravetools::FacePointer fp;
    fi=m.face.begin();
    for (i=0; i < m.fn; i++) {
      fp=&(*fi);
      if( ! fp->IsD() ) {
        itout(0,i) = indices[fp->cV(0)]+1;
        itout(1,i) = indices[fp->cV(1)]+1;
        itout(2,i) = indices[fp->cV(2)]+1;
      }
      ++fi;
    }
    return Rcpp::List::create(Rcpp::Named("vb") = vb,
                              Rcpp::Named("normals") = normals,
                              Rcpp::Named("it") = itout
    );

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
    return Rcpp::wrap(1);
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgUniformResample(
    const SEXP& vb_, const SEXP& it_, const float& voxelSize, const float& offsetThr,
    const bool& discretizeFlag, const bool& multiSampleFlag, const bool& absDistFlag,
    const bool& mergeCloseVert, const bool& silent) {
  try {
    ravetools::MyMesh m, baseMesh, offsetMesh;
    ravetools::IOMesh< ravetools::MyMesh >::vcgReadR( baseMesh , vb_ , it_ );
    if ( baseMesh.fn == 0 ) {
      Rcpp::stop( "This filter requires a mesh with some faces, it does not work on point cloud");
    }
    vcg::tri::UpdateBounding< ravetools::MyMesh >::Box(baseMesh);
    baseMesh.face.EnableNormal();
    vcg::Point3i volumeDim;
    vcg::Box3f volumeBox = baseMesh.bbox;
    volumeBox.Offset( volumeBox.Diag()/10.0f + offsetThr );

    BestDim(volumeBox , voxelSize, volumeDim );

    Rcpp::checkUserInterrupt();

    if (!silent) {
      Rprintf("Resampling mesh using a volume of %i x %i x %i\n", volumeDim[0], volumeDim[1], volumeDim[2] );
      Rprintf("  VoxelSize is %f, offset is %f\n", voxelSize, offsetThr);
      Rprintf("  Mesh Box is %f %f %f\n", baseMesh.bbox.DimX(),
              baseMesh.bbox.DimY(), baseMesh.bbox.DimZ() );
    }
    vcg::tri::Resampler< ravetools::MyMesh, ravetools::MyMesh >::Resample(
        baseMesh, offsetMesh, volumeBox, volumeDim, voxelSize * 3.5f,
        offsetThr, discretizeFlag, multiSampleFlag, absDistFlag
    );
    Rcpp::checkUserInterrupt();
    if ( mergeCloseVert ) {
      float mergeThr = offsetMesh.bbox.Diag() / 10000.0f;
      int total = vcg::tri::Clean< ravetools::MyMesh >::MergeCloseVertex( offsetMesh , mergeThr );
      if ( !silent ) {
        Rprintf("Successfully merged %d vertices with a distance lower than %f\n", total, mergeThr);
      }
    }
    vcg::tri::Allocator< ravetools::MyMesh >::CompactVertexVector( offsetMesh );
    vcg::tri::Allocator< ravetools::MyMesh >::CompactFaceVector( offsetMesh );
    vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertexAngleWeighted( offsetMesh );
    vcg::tri::UpdateNormal< ravetools::MyMesh >::NormalizePerVertex( offsetMesh );
    Rcpp::NumericMatrix vbout(3, offsetMesh.vn), normals(3, offsetMesh.vn);
    Rcpp::IntegerMatrix itout(3, offsetMesh.fn);
    vcg::SimpleTempData< ravetools::MyMesh::VertContainer, int > indiceout( offsetMesh.vert );
    ravetools::VertexIterator vi;
    vi = offsetMesh.vert.begin();
    for ( int i = 0 ; i < offsetMesh.vn ; i++ ) {
      indiceout[vi] = i;
      vbout(0,i) = (*vi).P()[0];
      vbout(1,i) = (*vi).P()[1];
      vbout(2,i) = (*vi).P()[2];
      normals(0,i) = (*vi).N()[0];
      normals(1,i) = (*vi).N()[1];
      normals(2,i) = (*vi).N()[2];
      ++vi;
    }
    Rcpp::checkUserInterrupt();
    ravetools::FaceIterator fi = offsetMesh.face.begin();
    for ( int i = 0; i < offsetMesh.fn ; i++, fi++ ) {
      itout(0, i) = indiceout[ fi->cV(0) ] + 1;
      itout(1, i) = indiceout[ fi->cV(1) ] + 1;
      itout(2, i) = indiceout[ fi->cV(2) ] + 1;
    }

    return Rcpp::List::create(Rcpp::Named("vb") = vbout,
                              Rcpp::Named("it") = itout,
                              Rcpp::Named("normals")=normals);

    return Rcpp::wrap(0);
  } catch ( std::exception& e ) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgUpdateNormals(SEXP vb_, SEXP it_, const int & select,
                      const Rcpp::IntegerVector & pointcloud, const bool & silent)
{
  try {
    ravetools::MyMesh m;
    ravetools::VertexIterator vi;

    // allocate mesh and fill it
    int check = ravetools::IOMesh< ravetools::MyMesh >::vcgReadR(m,vb_,it_);
    Rcpp::NumericMatrix normals(3, m.vn);
    if (check < 0) {
      Rcpp::stop("mesh has no faces and/or no vertices");
    } else if (check == 1) {
      if ( !silent ) {
        Rprintf("%s\n", "Info: mesh has no faces normals for point clouds are computed");
      }
      vcg::tri::PointCloudNormal< ravetools::MyMesh >::Param p;
      p.fittingAdjNum = pointcloud[0];
      p.smoothingIterNum = pointcloud[1];
      p.viewPoint = vcg::Point3f(0,0,0);
      p.useViewPoint = false;
      vcg::tri::PointCloudNormal< ravetools::MyMesh >::Compute(m,p);
    }  else {
      // update normals
      if (select == 0) {
        vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertex(m);
      } else {
        vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertexAngleWeighted(m);
      }
      vcg::tri::UpdateNormal< ravetools::MyMesh >::NormalizePerVertex(m);

      //write back
    }
    vi = m.vert.begin();
    vcg::SimpleTempData< ravetools::MyMesh::VertContainer , int > indiceout(m.vert);
    for ( int i = 0 ; i < m.vn ; i++) {
      if( ! vi->IsD() )	{
        normals(0,i) = (*vi).N()[0];
        normals(1,i) = (*vi).N()[1];
        normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }

    return Rcpp::wrap(normals);

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}

// [[Rcpp::export]]
SEXP vcgEdgeSubdivision(SEXP vb_, SEXP it_)
{
  try {
    ravetools::MyMesh m;

    // allocate mesh and fill it
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m,vb_,it_);

    // 1) build the Face -> Edge table
    std::vector<ravetools::MyPEdge> edgeVec;
    vcg::tri::UpdateTopology<ravetools::MyMesh>::FillUniqueEdgeVector(m, edgeVec, /*includeFauxEdge=*/false);

    // 2) remember original sizes
    const size_t origVn = m.vn;
    const size_t origEn = edgeVec.size();
    const size_t origFn = m.fn;

    // 3) allocate one new vertex PER EDGE
    vcg::tri::Allocator<ravetools::MyMesh>::AddVertices(m, origEn);

    // 4) compute each midpoint and build a map (edge->new-vertex)
    std::unordered_map<
      ravetools::EdgePair,
      ravetools::VertexPointer,
      ravetools::EdgeHash,
      ravetools::EdgeEqual> midMap;
    midMap.reserve(origEn);

    for(int i=0; i<origEn; ++i){
      ravetools::MyPEdge &pe = edgeVec[i];
      ravetools::FacePointer f = pe.f;               // face pointer
      int zi = pe.z;               // local index 0-2
      ravetools::VertexPointer v0 = f->V(zi);
      ravetools::VertexPointer v1 = f->V((zi+1)%3);
      // new midpoint vertex is at m.vert[origVn + i]
      ravetools::VertexPointer mv = &m.vert[ origVn + i ];
      mv->P() = (v0->P() + v1->P()) * ravetools::ScalarType(0.5);
      // AddVertices leaves the normal uninitialised; both parents have a
      // defined normal (see IOMesh::vcgReadR), so averaging is safe here
      mv->N() = (v0->N() + v1->N()) * ravetools::ScalarType(0.5);
      // canonical key ordering
      ravetools::EdgePair key = v0 < v1 ? ravetools::EdgePair(v0,v1) : ravetools::EdgePair(v1,v0);
      midMap[key] = mv;
    }

    // 5) now split each original face into four
    //    allocate 3 extra faces per original
    vcg::tri::Allocator<ravetools::MyMesh>::AddFaces(m, origFn*3);

    // 6) helper to get the midpoint-vertex for an edge
    for(int i=0;i<origFn;++i){
      ravetools::FacePointer f0 = &m.face[i];
      ravetools::VertexPointer A = f0->V(0), B = f0->V(1), C = f0->V(2);

      ravetools::EdgePair kAB = A < B ? ravetools::EdgePair(A,B) : ravetools::EdgePair(B,A);
      ravetools::EdgePair kBC = B < C ? ravetools::EdgePair(B,C) : ravetools::EdgePair(C,B);
      ravetools::EdgePair kCA = C < A ? ravetools::EdgePair(C,A) : ravetools::EdgePair(A,C);

      ravetools::VertexPointer MAB = midMap[kAB];
      ravetools::VertexPointer MBC = midMap[kBC];
      ravetools::VertexPointer MCA = midMap[kCA];

      // overwrite original face
      f0->V(0)=A;   f0->V(1)=MAB; f0->V(2)=MCA;

      // f1
      ravetools::FacePointer f1 = &m.face[ origFn + i ];
      f1->V(0)=MAB; f1->V(1)=B;   f1->V(2)=MBC;

      // f2
      ravetools::FacePointer f2 = &m.face[ origFn*2 + i ];
      f2->V(0)=MCA; f2->V(1)=MBC; f2->V(2)=C;

      // f3 (center)
      ravetools::FacePointer f3 = &m.face[ origFn*3 + i ];
      f3->V(0)=MAB; f3->V(1)=MBC; f3->V(2)=MCA;
    }

    Rcpp::List out = ravetools::IOMesh<ravetools::MyMesh>::vcgToR(m, false);
    return out;

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}

// [[Rcpp::export]]
SEXP vcgEdgeLengthSubdivision(SEXP vb_, SEXP it_, double maxEdgeLen, int maxIter)
{
  try {
    ravetools::MyMesh m;
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);

    m.face.EnableFFAdjacency();

    vcg::tri::MidPoint<ravetools::MyMesh> midFunctor(&m);
    vcg::tri::EdgeLen<ravetools::MyMesh, ravetools::ScalarType> edgePred(
        static_cast<ravetools::ScalarType>(maxEdgeLen));

    for (int iter = 0; iter < maxIter; iter++) {
      vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);
      bool refined = vcg::tri::RefineE<
          ravetools::MyMesh,
          vcg::tri::MidPoint<ravetools::MyMesh>,
          vcg::tri::EdgeLen<ravetools::MyMesh, ravetools::ScalarType>>(
          m, midFunctor, edgePred);
      if (!refined) break;
    }

    vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(m);
    return ravetools::IOMesh<ravetools::MyMesh>::vcgToR(m, false);

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}

// [[Rcpp::export]]
SEXP vcgVolume( SEXP mesh_ )
{
  try {
    ravetools::MyMesh m;
    ravetools::IOMesh<ravetools::MyMesh>::mesh3d2vcg(m, mesh_);
    bool Watertight, Oriented = false;
    int VManifold, FManifold;
    float Volume = 0;
    // int numholes, BEdges = 0;
    //check manifoldness
    m.vert.EnableVFAdjacency();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    m.face.EnableNormal();
    vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);
    VManifold = vcg::tri::Clean<ravetools::MyMesh>::CountNonManifoldVertexFF(m);
    FManifold = vcg::tri::Clean<ravetools::MyMesh>::CountNonManifoldEdgeFF(m);

    if ((VManifold>0) || (FManifold>0)) {
      throw std::runtime_error(
        (
            "Mesh is not manifold\n  Non-manifold vertices: " +
              std::to_string(VManifold) +"\n" +
              "  Non-manifold edges: " +
              std::to_string(FManifold) +"\n"
        ).c_str()
      );
    }


    Watertight = vcg::tri::Clean<ravetools::MyMesh>::IsWaterTight(m);
    Oriented = vcg::tri::Clean<ravetools::MyMesh>::IsCoherentlyOrientedMesh(m);
    vcg::tri::Inertia<ravetools::MyMesh> mm(m);
    mm.Compute(m);
    Volume = mm.Mass();

    // the sign of the volume depend on the mesh orientation
    if (Volume < 0.0)
      Volume = -Volume;
    if (!Watertight)
      ::Rf_warning("Mesh is not watertight! USE RESULT WITH CARE!\n");
    if (!Oriented)
      ::Rf_warning("Mesh is not coherently oriented! USE RESULT WITH CARE!\n");

    return Rcpp::wrap(Volume);

  } catch (std::exception& e) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}


// [[Rcpp::export]]
SEXP vcgSphere(const int& subdiv, bool normals) {
  try {
    ravetools::MyMesh m;
    Sphere(m,subdiv);
    if (normals)
      vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexNormalized(m);
    Rcpp::List out = ravetools::IOMesh<ravetools::MyMesh>::vcgToR(m,normals);
    return out;
  } catch (std::exception& e) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgDijkstra(SEXP vb_, SEXP it_, const Rcpp::IntegerVector & source, const double & maxdist_) {
  try {

    ravetools::ScalarType maxdist = std::numeric_limits<ravetools::ScalarType>::max();
    if( maxdist_ != NA_REAL && maxdist_ > 0.0 ) {
      maxdist = maxdist_;
    }

    // Declare Mesh and helper variables
    ravetools::MyMesh m;
    ravetools::VertexIterator vi;

    // Allocate mesh and fill it
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.vert.EnableQuality();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    vcg::tri::UpdateTopology<ravetools::MyMesh>::VertexFace(m);

    // Create int vertex indices to return to R.
    vcg::SimpleTempData<ravetools::MyMesh::VertContainer,int> indices(m.vert);
    vi = m.vert.begin();
    for (int i=0; i < m.vn; i++) {
      indices[vi] = i;
      ++vi;
    }

    // Prepare seed vector with source vertex
    std::vector<ravetools::MyVertex*> seedVec;
    for ( int i = 0; i < source.length(); i++ ) {
      vi = m.vert.begin() + source[i];
      seedVec.push_back( &*vi );
    }

    std::vector<ravetools::MyVertex*> inInterval;
    ravetools::MyMesh::PerVertexAttributeHandle<ravetools::MyMesh::VertexPointer> sourcesHandle;
    sourcesHandle = vcg::tri::Allocator<ravetools::MyMesh>::AddPerVertexAttribute<ravetools::MyMesh::VertexPointer> (m, "sources");
    ravetools::MyMesh::PerVertexAttributeHandle<ravetools::MyMesh::VertexPointer> parentHandle;
    parentHandle = vcg::tri::Allocator<ravetools::MyMesh>::AddPerVertexAttribute<ravetools::MyMesh::VertexPointer> (m, "parent");

    // Compute pseudo-geodesic distance by summing dists along shortest path in graph.
    vcg::tri::EuclideanDistance<ravetools::MyMesh> ed;
    vcg::tri::Geodesic<ravetools::MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed, maxdist, &inInterval, &sourcesHandle, &parentHandle);
    std::vector<double> geodist;
    std::vector<int> parentNode;

    ravetools::MyMesh::VertexPointer parent;
    // ravetools::MyMesh::VertexPointer source;
    int parentidx = 0;

    vi=m.vert.begin();
    for ( int i = 0 ; i < m.vn ; i++, vi++ ) {
      parent = parentHandle[ i ];
      // source = sourcesHandle[ i ];
      if( parent == NULL ) {
        // Rcpp::Rcout << i << " -> NA";
        parentNode.push_back( NA_INTEGER );
        geodist.push_back( NA_REAL );
      } else {
        parentidx = indices[parent];
        // Rcpp::Rcout << i << " -> " << parentidx;
        if( parentidx == i ) {
          // source node
          parentNode.push_back( NA_INTEGER );
        } else {
          parentNode.push_back( parentidx );
        }
        geodist.push_back( (double) ( vi->Q() ) );
      }
      // if( source != NULL ) {
      //   Rcpp::Rcout << "  src: " << indices[source];
      // }
      // Rcpp::Rcout << "\n";
    }

    // clean up
    vcg::tri::Allocator<ravetools::MyMesh>::DeletePerVertexAttribute<ravetools::MyMesh::VertexPointer> (m, sourcesHandle);
    vcg::tri::Allocator<ravetools::MyMesh>::DeletePerVertexAttribute<ravetools::MyMesh::VertexPointer> (m, parentHandle);

    // parents are 0-indexed
    Rcpp::List L = Rcpp::List::create(Rcpp::Named("parent") = parentNode , Rcpp::Named("geodist") = geodist);
    return L;
  } catch (std::exception& e) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}


// [[Rcpp::export]]
SEXP vcgRaycaster(
    SEXP vb_ , SEXP it_,
    const Rcpp::NumericVector & rayOrigin, // 3 x n matrix
    const Rcpp::NumericVector & rayDirection,
    const float & maxDistance,
    const bool & bothSides,
    const int & threads = 1)
{
  try {
    ravetools::MyMesh m;
    int check = ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);
    if (check != 0) {
      throw std::runtime_error("Mesh has no faces or vertices. Unable to perform raycaster");
    }

    ravetools::ScalarType x,y,z;
    ravetools::MyMesh rays;

    // Leave R to check if nRays > 0
    unsigned int nRays = rayOrigin.length() / 3;
    Rcpp::NumericVector castDistance(nRays);
    Rcpp::IntegerVector hitFlag(nRays);

    // rayOrigin and rayDirection must be identical 3xn
    // Leave the checks in R wrapper
    Rcpp::NumericVector intersectPoints(nRays * 3);
    Rcpp::NumericVector intersectNormals(nRays * 3);
    Rcpp::IntegerVector intersectIndex(nRays);

    //Allocate target
    std::vector<ravetools::MyMesh::VertexPointer> ivp;
    vcg::tri::Allocator<ravetools::MyMesh>::AddVertices(rays, nRays);
    vcg::Point3f normtmp;

    // Copy the rayOrigin and rayDirection
    ravetools::MyMesh::VertexIterator vi = rays.vert.begin();
    for (unsigned int i=0; i < nRays; i++, vi++) {
      x = rayOrigin[ i * 3 ];
      y = rayOrigin[ i * 3 + 1 ];
      z = rayOrigin[ i * 3 + 2 ];
      (*vi).P() = ravetools::MyMesh::CoordType(x, y, z);
      x = rayDirection[ i * 3 ];
      y = rayDirection[ i * 3 + 1 ];
      z = rayDirection[ i * 3 + 2 ];
      // Rcpp::Rcout << x << " " << y << " " << z << "\n";
      (*vi).N() = ravetools::MyMesh::CoordType(x, y, z);
    }

    // bounding box to calculate max cast distance
    m.face.EnableNormal();
    vcg::tri::UpdateBounding<ravetools::MyMesh>::Box(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::PerFaceNormalized(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(m);

    vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(rays);

    vcg::tri::FaceTmark<ravetools::MyMesh> mf;
    mf.SetMesh( &m );
    vcg::RayTriangleIntersectionFunctor<true> FintFunct;
    vcg::GridStaticPtr<ravetools::MyMesh::FaceType, ravetools::MyMesh::ScalarType> gridSearcher;
    gridSearcher.Set(m.face.begin(), m.face.end());

#pragma omp parallel for firstprivate(maxDistance,gridSearcher,mf) schedule(static) num_threads(threads)
{
    for ( int i = 0; i < rays.vn ; i++ ) {
      float t0 = 0.0f, t1 = 0.0f;
      int faceIndex = -1;
      vcg::Ray3f ray;
      vcg::Point3f orig = rays.vert[i].P();
      // vcg::Point3f orig0 = orig;
      vcg::Point3f dir = rays.vert[i].N();
      vcg::Point3f intersection = ravetools::MyMesh::CoordType(0, 0, 0);
      ravetools::MyFace* facePtr0;
      ravetools::MyFace* facePtr1;

      /**
       *  Set ray origin to be slightly "behind" the `orig`
       *  This is because if orig coincide with the underlying intersection,
       *  FintFunct will not be able to identify the intersection
       *  Example:
       sphere <- ravetools::vcg_sphere()
       box <- Rvcg::vcgBox(sphere)
       box$vb[1:3,] <- box$vb[1:3,] + c(1,1,1) - 1e-6
       mesh <- box
       vcgRaycaster(vb_ = mesh$vb, it_ = mesh$it - 1L, rayOrigin = matrix(c(0,0,0), ncol = 1), rayDirection = matrix(c(1,1,1), ncol = 1), maxDistance = 1e14, bothSides = FALSE)
       */
      ray.SetOrigin(orig - 1e-6f * dir);
      ray.SetDirection(dir);

      // raycaster
      facePtr0 = GridDoRay(gridSearcher, FintFunct, mf, ray, maxDistance, t0);
      if ( bothSides ) {
        ray.SetOrigin(orig + 1e-6f * dir);
        // cast the ray backwards
        ray.SetDirection(-dir);
        facePtr1 = GridDoRay(gridSearcher, FintFunct, mf, ray, maxDistance, t1);
        if( facePtr1 && ( !facePtr0 || t1 < t0 ) ) {
          facePtr0 = facePtr1;
          t0 = -t1;
        }
      }

      if( facePtr0 ) {
        // pay off the debt
        if( t0 > 0.0f ) {
          t0 -= 1e-6f;
        } else {
          t0 += 1e-6f;
        }

        intersection = rays.vert[i].P()+dir * t0;
        castDistance[ i ] = t0;
        hitFlag[ i ] = 1;

        faceIndex = vcg::tri::Index(m, facePtr0);

        // face normal
        const ravetools::MyFace face = m.face[faceIndex];
        ravetools::MyMesh::CoordType faceNormal = (face.V(0)->N() + face.V(1)->N() + face.V(2)->N()).normalized();

        intersectPoints[i * 3] = intersection[0];
        intersectPoints[i * 3 + 1] = intersection[1];
        intersectPoints[i * 3 + 2] = intersection[2];
        intersectNormals[i * 3] = faceNormal[0];
        intersectNormals[i * 3 + 1] = faceNormal[1];
        intersectNormals[i * 3 + 2] = faceNormal[2];
        intersectIndex[i] = faceIndex;

      } else {
        // No intersection
        castDistance[ i ] = NA_REAL;
        hitFlag[ i ] = 0;

        intersectPoints[i * 3] = NA_REAL;
        intersectPoints[i * 3 + 1] = NA_REAL;
        intersectPoints[i * 3 + 2] = NA_REAL;
        intersectNormals[i * 3] = NA_REAL;
        intersectNormals[i * 3 + 1] = NA_REAL;
        intersectNormals[i * 3 + 2] = NA_REAL;
        intersectIndex[i] = NA_INTEGER;
      }
    }
}
    return Rcpp::List::create(Rcpp::Named("intersectPoints") = intersectPoints,
                              Rcpp::Named("intersectNormals") = intersectNormals,
                              Rcpp::Named("intersectIndex") = intersectIndex,
                              Rcpp::Named("hitFlag") = hitFlag,
                              Rcpp::Named("castDistance") = castDistance
    );
  } catch (std::exception& e) {
    Rcpp::stop( e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue;
}


// [[Rcpp::export]]
SEXP vcgKDTreeSearch(
    SEXP target_, SEXP query_,
    unsigned int k,
    unsigned int nPointsPerCell = 16,
    unsigned int maxDepth = 64
)
{
  try {

    ravetools::MyPointCloud target, query;
    ravetools::IOMesh<ravetools::MyPointCloud>::vcgReadR(target, target_);
    ravetools::IOMesh<ravetools::MyPointCloud>::vcgReadR(query, query_);

    // List out = Rvcg::KDtree< PcMesh, PcMesh >::KDtreeIO(target, query, k,nofP, mDepth,threads);
    // typedef std::pair<float,int> mypair;
    Rcpp::IntegerMatrix index(query.vn, k);
    Rcpp::NumericMatrix distance(query.vn, k);
    std::fill(index.begin(), index.end(), -1);

    vcg::VertexConstDataWrapper<ravetools::MyPointCloud> targetWrapper(target);
    vcg::KdTree<float> tree(targetWrapper, nPointsPerCell, maxDepth);

    //tree.setMaxNofNeighbors(k);
    vcg::KdTree<float>::PriorityQueue queue;

    std::vector< std::pair<float,int> > sortPairs;

    for (int i = 0; i < query.vn; i++) {
      tree.doQueryK(query.vert[i].cP(), k, queue);
      int neighbors = queue.getNofElements();

      sortPairs.clear();

      for (int j = 0; j < neighbors; j++) {
        int neightId = queue.getIndex(j);
        float dist = Distance(query.vert[i].cP(), target.vert[neightId].cP());
        sortPairs.push_back( std::pair<float,int>(dist, neightId) );
      }

      std::sort(sortPairs.begin(), sortPairs.end());
      for (int j = 0; j < neighbors; j++){
        index(i, j) = sortPairs[j].second;
        distance(i, j) = sortPairs[j].first;
      }
    }

    return Rcpp::List::create(
      Rcpp::Named("index") = index,
      Rcpp::Named("distance") = distance
    );

  } catch (std::exception& e) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}

// [[Rcpp::export]]
SEXP vcgSubset(SEXP vb_ , SEXP it_, const Rcpp::LogicalVector selector_) {
  try {
    // allocate mesh and fill it
    // Load mesh from R data
    ravetools::MyMesh m;
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);

    // Mark vertices for deletion based on selector_
    size_t vertCount = m.vert.size();

    if (selector_.length() != vertCount) {
      Rcpp::stop("Inconsistent lengths");
    }

    if (selector_.size() != static_cast<int>(vertCount)) {
      Rcpp::stop("Selector length does not match number of vertices");
    }
    for (size_t i = 0; i < vertCount; ++i) {
      if (selector_[i]) {
        m.vert[i].SetD();  // mark vertex for deletion
      }
    }

    // Mark faces for deletion if any of their vertices was deleted
    for (auto &f : m.face) {
      if (f.IsD()) continue;  // skip already deleted faces
      for (int vi = 0; vi < 3; ++vi) {
        if (f.V(vi)->IsD()) {
          f.SetD();          // delete the face
          break;
        }
      }
    }

    vcg::tri::Allocator<ravetools::MyMesh>::CompactVertexVector(m);
    vcg::tri::Allocator<ravetools::MyMesh>::CompactFaceVector(m);

    // Return the subset mesh back to R
    // vcgToR will take care of re-ordering
    Rcpp::List out = ravetools::IOMesh<ravetools::MyMesh>::vcgToR(m, false);
    return out;

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall

}

// [[Rcpp::export]]
Rcpp::List vcgCountEdgeDefects(SEXP vb_, SEXP it_) {
  try {
    ravetools::MyMesh m;
    int check = ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);
    if (check < 0) {
      Rcpp::stop("vcgCountEdgeDefects: mesh has no faces and/or no vertices");
    }

    int boundary = 0, nonmanifold = 0;
    countEdgeDefects(m, boundary, nonmanifold);

    return Rcpp::List::create(
      Rcpp::Named("boundary_edges")    = boundary,
      Rcpp::Named("nonmanifold_edges") = nonmanifold,
      Rcpp::Named("is_closed_manifold") = (boundary == 0 && nonmanifold == 0)
    );

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}

// [[Rcpp::export]]
double vcgAverageEdgeLength(SEXP vb_, SEXP it_) {
  try {
    ravetools::MyMesh m;
    int check = ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);
    if (check < 0) {
      Rcpp::stop("vcgAverageEdgeLength: mesh has no faces and/or no vertices");
    }

    return averageEdgeLength(m);

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return NA_REAL; // -Wall
}

// [[Rcpp::export]]
double vcgMaxEdgeLength(SEXP vb_, SEXP it_) {
  try {
    ravetools::MyMesh m;
    int check = ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);
    if (check < 0) {
      Rcpp::stop("vcgMaxEdgeLength: mesh has no faces and/or no vertices");
    }
    return maxEdgeLength(m);
  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return NA_REAL; // -Wall
}

// [[Rcpp::export]]
Rcpp::IntegerVector vcgMeshPatchFaces(SEXP vb_, SEXP it_,
                                       Rcpp::IntegerVector boundary_seq,
                                       int seed_face)
{
  try {
    ravetools::MyMesh m;
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);
    int nv = (int)m.VN();
    int nf = (int)m.FN();

    m.face.EnableFFAdjacency();
    vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);

    // Build boundary edge set (canonical key: min*nv + max)
    std::unordered_set<int64_t> cut_edges;
    cut_edges.reserve(boundary_seq.size());
    int nb = (int)boundary_seq.size();
    for (int k = 0; k < nb; k++) {
      int a = boundary_seq[k], b = boundary_seq[(k + 1) % nb];
      cut_edges.insert((int64_t)std::min(a, b) * nv + std::max(a, b));
    }

    // BFS: returns 1-based face indices reachable from seed_fi
    auto bfs = [&](int seed_fi) -> std::vector<int> {
      std::vector<bool> vis(nf, false);
      std::vector<int>  q(nf);
      int head = 0, tail = 0;
      vis[seed_fi] = true;
      q[tail++] = seed_fi;
      while (head < tail) {
        int fi = q[head++];
        ravetools::FacePointer fp = &m.face[fi];
        for (int e = 0; e < 3; e++) {
          if (fp->FFp(e) == fp) continue; // mesh boundary edge
          int va  = (int)vcg::tri::Index(m, fp->V(e));
          int vbi = (int)vcg::tri::Index(m, fp->V((e + 1) % 3));
          if (cut_edges.count((int64_t)std::min(va, vbi) * nv + std::max(va, vbi)))
            continue; // cut edge
          int nb_fi = (int)vcg::tri::Index(m, fp->FFp(e));
          if (!vis[nb_fi]) { vis[nb_fi] = true; q[tail++] = nb_fi; }
        }
      }
      std::vector<int> result;
      result.reserve(tail);
      for (int i = 0; i < nf; i++)
        if (vis[i] && !m.face[i].IsD()) result.push_back(i + 1); // 1-based
      return result;
    };

    // Determine seed face
    int seed_fi  = seed_face;
    int other_fi = -1;

    if (seed_face < 0) {
      // Find two faces adjacent to first boundary edge
      int v0 = boundary_seq[0], v1 = boundary_seq[1];
      int64_t key = (int64_t)std::min(v0, v1) * nv + std::max(v0, v1);
      for (auto fi = m.face.begin(); fi != m.face.end() && seed_fi < 0; ++fi) {
        if (fi->IsD()) continue;
        for (int e = 0; e < 3; e++) {
          int va  = (int)vcg::tri::Index(m, fi->V(e));
          int vbi = (int)vcg::tri::Index(m, fi->V((e + 1) % 3));
          if ((int64_t)std::min(va, vbi) * nv + std::max(va, vbi) == key) {
            seed_fi = (int)vcg::tri::Index(m, &*fi);
            ravetools::FacePointer nbfp = fi->FFp(e);
            if (nbfp != &*fi)
              other_fi = (int)vcg::tri::Index(m, nbfp);
            break;
          }
        }
      }
      if (seed_fi < 0)
        Rcpp::stop("vcgMeshPatchFaces: could not find a face adjacent to the first boundary edge.");
      if (other_fi < 0)
        Rcpp::stop("vcgMeshPatchFaces: first boundary edge is on the mesh boundary; please specify seed_vertex explicitly.");
    }

    std::vector<int> patch = bfs(seed_fi);

    // Non-dividing: loop does not partition the mesh, return all faces
    if ((int)patch.size() == nf) return Rcpp::wrap(patch);

    // Auto-select smaller side
    if (seed_face < 0 && (int)patch.size() > nf / 2) {
      patch = bfs(other_fi);
    }

    return Rcpp::wrap(patch);

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return Rcpp::IntegerVector(); // -Wall
}

// ---------------------------------------------------------------------------
// Detect and repair common defects in triangular surface meshes so they
// become closed, manifold, genus-0 surfaces, a hard precondition for
// algorithms such as mris_inflate (see mrisCheckClosedManifold in mrisCommon.cpp).
//
// Typical source of these defects: isosurfaces extracted from volumes
// (e.g. vcg_isosurface) via marching-cubes-style algorithms, which can leave
// behind small "cracks", isolated boundary-edge loops bounding tiny holes
// that are not caused by duplicated/coincident vertices (so simple vertex
// welding does not close them) but are genuine small gaps in the tessellation.
//
// Repair strategy (all via VCG's existing, well-tested algorithms):
//   1. Remove degenerate / duplicate faces.
//   2. Weld near-coincident vertices (closes cracks that *are* caused by
//      duplicated vertices).
//   3. Ear-cut fill any remaining small boundary loops ("isolated edges").
//   4. Remove unreferenced vertices and compact.
//   5. Re-orient faces coherently (consistent winding order across the
//      surface) and, if the result is a single watertight component, flip
//      normals to point outward.
//
// Parameters:
//   merge_tolerance  distance (in mesh units) below which vertices are welded
//                    together; non-positive => auto (1e-4 * average edge length)
//   max_hole_size    maximum number of boundary edges of a hole that will be
//                    triangulated (ear-cutting fill); holes larger than this
//                    are left untouched (and will be reported as remaining
//                    boundary edges)
//   verbose          print a short before/after diagnostic report
//
// Returns a list with the repaired mesh (vb, it, normals) plus an `info`
// list describing what was found/changed.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List vcgFixDefects(
    SEXP   vb_,
    SEXP   it_,
    double merge_tolerance = -1.0,
    int    max_hole_size   = 100,
    bool   verbose         = false
) {
    try {
        ravetools::MyMesh m;
        int check = ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m, vb_, it_);
        if (check < 0) {
            Rcpp::stop("vcgFixDefects: mesh has no faces and/or no vertices");
        }

        m.vert.EnableVFAdjacency();
        m.face.EnableFFAdjacency();
        m.face.EnableVFAdjacency();
        m.face.EnableNormal();

        vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);
        int boundary_before = 0, nonmanifold_before = 0;
        countEdgeDefects(m, boundary_before, nonmanifold_before);

        if (verbose) {
            Rprintf("vcgFixDefects: input nv=%d nf=%d, boundary edges=%d, non-manifold edges=%d\n",
                    m.vn, m.fn, boundary_before, nonmanifold_before);
        }

        // 1. Remove degenerate / duplicate faces (zero-area, repeated triangles)
        vcg::tri::Clean<ravetools::MyMesh>::RemoveDegenerateFace(m);
        Rcpp::checkUserInterrupt();

        vcg::tri::Clean<ravetools::MyMesh>::RemoveDuplicateFace(m);
        Rcpp::checkUserInterrupt();

        if (verbose) {
            Rprintf("vcgFixDefects: [1] removed degenerate/duplicate faces -> nv=%d nf=%d\n",
                    m.vn, m.fn);
        }

        // 2. Weld near-coincident vertices, closes cracks that arise from
        //    duplicated vertices at (near-)identical positions.
        double tol = merge_tolerance;
        if (tol <= 0.0) {
            double avgLen = averageEdgeLength(m);
            tol = avgLen * 1e-4;
        }

        int merged = 0;
        if (tol > 0.0) {
            merged = vcg::tri::Clean<ravetools::MyMesh>::MergeCloseVertex(
                m, (ravetools::ScalarType)tol);
        }

        Rcpp::checkUserInterrupt();
        if (verbose) {
            Rprintf("vcgFixDefects: [2] merged %d close vertices -> nv=%d nf=%d\n",
                    merged, m.vn, m.fn);
        }

        // 3. Ear-cut fill any remaining small boundary loops ("isolated
        //    edges" / small holes that are genuine gaps, not duplicate-vertex
        //    cracks, and therefore are not closed by welding).
        // Rebuild FF adjacency *before* compacting. Steps 1-2 above mark faces
        // deleted (RemoveDegenerateFace / RemoveDuplicateFace / MergeCloseVertex)
        // without detaching the FF links of the surviving faces that point at
        // them. CompactFaceVector then remaps every FF link through
        // PointerUpdater::remap, whose entries for deleted faces are left at
        // size_t(-1), so it evaluates `fbase + size_t(-1)` -- pointer overflow
        // (clang-UBSAN: allocate.h:1330) that stores a wild adjacency pointer
        // one element before the face array. Rebuilding first guarantees every
        // FF link refers to a surviving face; the compactor then remaps them
        // correctly and topology stays valid for the hole fill below.
        vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);
        vcg::tri::Allocator<ravetools::MyMesh>::CompactEveryVector(m);
        Rcpp::checkUserInterrupt();

        // MinimumWeightEar inspects neighboring faces' and vertices' normals
        // (FFlip()->cN(), e0.v->N()) while picking the best ear to cut, so
        // both must be allocated *and* populated before EarCuttingFill runs.
        vcg::tri::UpdateNormal<ravetools::MyMesh>::PerFaceNormalized(m);
        vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexNormalized(m);
        Rcpp::checkUserInterrupt();

        if (verbose) {
          Rprintf("vcgFixDefects: [3a] topology/normals ready, starting hole fill (max_hole_size=%d)\n", max_hole_size);
        }
        int holes_filled = vcg::tri::Hole<ravetools::MyMesh>::template EarCuttingFill<
            vcg::tri::MinimumWeightEar<ravetools::MyMesh> >(m, max_hole_size, false);
        if (verbose) {
          Rprintf("vcgFixDefects: [3b] filled %d hole(s) -> nv=%d nf=%d\n",
                  holes_filled, m.vn, m.fn);
        }

        // 4. Remove unreferenced vertices, compact containers
        vcg::tri::Clean<ravetools::MyMesh>::RemoveUnreferencedVertex(m);
        vcg::tri::Allocator<ravetools::MyMesh>::CompactEveryVector(m);
        Rcpp::checkUserInterrupt();
        if (verbose) {
          Rprintf("vcgFixDefects: [4] removed unreferenced vertices -> nv=%d nf=%d\n", m.vn, m.fn);
        }

        // 5. Fix face winding order: orient all faces coherently
        vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);
        if (verbose) {
          Rprintf("vcgFixDefects: [5a] topology rebuilt, orienting coherently\n");
        }
        bool is_oriented = false, is_orientable = false;
        vcg::tri::Clean<ravetools::MyMesh>::OrientCoherentlyMesh(m, is_oriented, is_orientable);
        Rcpp::checkUserInterrupt();
        if (verbose) {
          Rprintf("vcgFixDefects: [5b] oriented=%d orientable=%d\n",
                  (int)is_oriented, (int)is_orientable);
        }

        // If the mesh is now a single watertight, coherently-oriented
        // component, make sure normals point outward (assumes watertight,
        // as documented by VCG's FlipNormalOutside).
        bool flipped = false;
        int boundary_after = 0, nonmanifold_after = 0;
        countEdgeDefects(m, boundary_after, nonmanifold_after);
        if (is_orientable && boundary_after == 0 && nonmanifold_after == 0) {
            vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);
            flipped = vcg::tri::Clean<ravetools::MyMesh>::FlipNormalOutside(m);
        }
        Rcpp::checkUserInterrupt();

        // Recompute normals for output
        vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexNormalizedPerFace(m);
        Rcpp::checkUserInterrupt();

        vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(m);
        Rcpp::checkUserInterrupt();

        if (verbose) {
            Rprintf("vcgFixDefects: removed/merged %d vertices (tol=%.6g), filled %d hole(s)\n",
                    merged, tol, holes_filled);
            Rprintf("vcgFixDefects: output nv=%d nf=%d, boundary edges=%d, non-manifold edges=%d, "
                    "oriented=%s, orientable=%s, normals_flipped_outward=%s\n",
                    m.vn, m.fn, boundary_after, nonmanifold_after,
                    is_oriented ? "yes" : "no", is_orientable ? "yes" : "no",
                    flipped ? "yes" : "no");
        }

        Rcpp::List out = ravetools::IOMesh<ravetools::MyMesh>::vcgToR(m, true);
        out["info"] = Rcpp::List::create(
            Rcpp::Named("boundary_edges_before")    = boundary_before,
            Rcpp::Named("boundary_edges_after")     = boundary_after,
            Rcpp::Named("nonmanifold_edges_before") = nonmanifold_before,
            Rcpp::Named("nonmanifold_edges_after")  = nonmanifold_after,
            Rcpp::Named("vertices_merged")          = merged,
            Rcpp::Named("merge_tolerance")          = tol,
            Rcpp::Named("holes_filled")             = holes_filled,
            Rcpp::Named("is_oriented")              = is_oriented,
            Rcpp::Named("is_orientable")            = is_orientable,
            Rcpp::Named("normals_flipped_outward")  = flipped,
            Rcpp::Named("is_closed_manifold")       = (boundary_after == 0 && nonmanifold_after == 0)
        );

        return out;

    } catch (std::exception& e) {
        Rcpp::stop(e.what());
    } catch (...) {
        Rcpp::stop("unknown exception");
    }
    return R_NilValue; // -Wall
}
