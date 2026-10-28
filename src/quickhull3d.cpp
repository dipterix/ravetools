// vcg_convex_hull(method = "quickhull3d"): runs the port in quickhull3d.h. It
// is a translation unit of its own because that header turns off clang's
// floating-point contraction for everything that follows it.

#include <Rcpp.h>
#include "quickhull3d.h"

// [[Rcpp::export]]
Rcpp::IntegerMatrix quickhull3dConvexHull(const Rcpp::NumericMatrix &vb,
                                          double tolerance)
{
  try {
    if (vb.nrow() != 3) {
      Rcpp::stop("vcg_convex_hull: points must be a matrix with 3 rows");
    }
    ravetools::quickhull3d::QuickHull hull(vb.begin(), vb.ncol(), tolerance);
    hull.build();
    std::vector<std::array<int, 3> > faces = hull.collectFaces();

    // 1-based indices of the input points
    Rcpp::IntegerMatrix out(3, faces.size());
    for (std::size_t k = 0; k < faces.size(); k++) {
      for (int j = 0; j < 3; j++) {
        out(j, k) = faces[k][j] + 1;
      }
    }
    return out;

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return Rcpp::IntegerMatrix(); // -Wall
}
