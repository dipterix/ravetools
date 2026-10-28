// Convex hull of a point cloud (vcg_convex_hull).
//
// vcgConvexHull wraps vcglib's Quickhull (vcg/complex/algorithms/convex_hull.h).
// The point and hull types are separate from ravetools::MyMesh: the hull
// needs face-face adjacency and face normals as plain components, and double
// precision keeps the visibility tests consistent on the many coplanar points
// of voxel-derived inputs, where float rounding flips them. The result is the
// hull faces as 1-based indices of the input points; R/vcg-convex-hull.R
// validates the input and builds the mesh.

#include <Rcpp.h>
#include "vcgCommon.h"
#include <vcg/complex/algorithms/convex_hull.h>

namespace {

class HullPointVertex;
struct HullPointUsedTypes : public vcg::UsedTypes<
  vcg::Use<HullPointVertex>::AsVertexType
> {};

// the normal is unused by the hull, but IOMesh::vcgReadR zeroes a Point3d
// normal on every vertex it reads (vcglib's placeholder normal is a Point3f)
class HullPointVertex : public vcg::Vertex<
  HullPointUsedTypes,
  vcg::vertex::Coord3d,
  vcg::vertex::Normal3d,
  vcg::vertex::BitFlags
> {};

class HullPointCloud : public vcg::tri::TriMesh<
  std::vector<HullPointVertex>
> {};

class HullVertex;
class HullFace;
struct HullUsedTypes : public vcg::UsedTypes<
  vcg::Use<HullVertex>::AsVertexType,
  vcg::Use<HullFace>::AsFaceType
> {};

class HullVertex : public vcg::Vertex<
  HullUsedTypes,
  vcg::vertex::Coord3d,
  vcg::vertex::Normal3d,
  vcg::vertex::BitFlags
> {};

class HullFace : public vcg::Face<
  HullUsedTypes,
  vcg::face::VertexRef,
  vcg::face::FFAdj,
  vcg::face::Normal3d,
  vcg::face::BitFlags
> {};

class HullMesh : public vcg::tri::TriMesh<
  std::vector<HullVertex>, std::vector<HullFace>
> {};

} // namespace

// [[Rcpp::export]]
Rcpp::IntegerMatrix vcgConvexHull(SEXP vb_, double tolerance)
{
  try {
    HullPointCloud points;
    // 1 (no faces) is the expected result for a point cloud
    if (ravetools::IOMesh<HullPointCloud>::vcgReadR(points, vb_) < 0) {
      Rcpp::stop("vcg_convex_hull: points must be a matrix with 3 rows");
    }

    HullMesh hull;
    if (!vcg::tri::ConvexHull<HullPointCloud, HullMesh>::ComputeConvexHull(
        points, hull, tolerance)) {
      Rcpp::stop("vcg_convex_hull: at least 4 points are needed");
    }

    // The vertices ComputeConvexHull discards stay in hull.vert, only flagged
    // deleted, so the faces are mapped to the input through the attribute it
    // records rather than through the hull's own vertex order
    HullMesh::PerVertexAttributeHandle<size_t> indexInput =
      vcg::tri::Allocator<HullMesh>::GetPerVertexAttribute<size_t>(
        hull, std::string("indexInput"));

    Rcpp::IntegerMatrix faces(3, hull.fn);
    int k = 0;
    for (HullMesh::FaceIterator fi = hull.face.begin(); fi != hull.face.end(); ++fi) {
      if (fi->IsD()) continue;
      for (int j = 0; j < 3; j++) {
        faces(j, k) = (int)indexInput[fi->V(j)] + 1;
      }
      k++;
    }
    return faces;

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return Rcpp::IntegerMatrix(); // -Wall
}
