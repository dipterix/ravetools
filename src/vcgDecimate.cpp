// Quadric edge-collapse decimation (vcg_decimate).
//
// Adapted from Stefan Schlager's `Rvcg` (src/RQEdecim.cpp, GPL >= 2), itself
// an adaption of the `tridecimator` example shipped with vcglib. The mesh type
// is separate from ravetools::MyMesh because the collapse needs a per-vertex
// quadric and plain (non-optional) vertex-face adjacency.

#include <Rcpp.h>
#include "vcgCommon.h"
#include <vcg/math/quadric.h>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

namespace {

class DecVertex;
class DecEdge;
class DecFace;

struct DecUsedTypes: public vcg::UsedTypes<
  vcg::Use<DecVertex>::AsVertexType,
  vcg::Use<DecEdge>::AsEdgeType,
  vcg::Use<DecFace>::AsFaceType
> {};

class DecVertex : public vcg::Vertex<
  DecUsedTypes,
  vcg::vertex::VFAdj,
  vcg::vertex::Coord3f,
  vcg::vertex::Normal3f,
  vcg::vertex::Mark,
  vcg::vertex::BitFlags
> {
public:
  vcg::math::Quadric<double> &Qd() { return q; }
private:
  vcg::math::Quadric<double> q;
};

class DecEdge : public vcg::Edge<DecUsedTypes> {};

class DecFace : public vcg::Face<
  DecUsedTypes,
  vcg::face::VFAdj,
  vcg::face::VertexRef,
  vcg::face::BitFlags
> {};

class DecMesh : public vcg::tri::TriMesh<
  std::vector<DecVertex>, std::vector<DecFace>
> {};

typedef vcg::tri::BasicVertexPair<DecVertex> DecVertexPair;

class DecCollapse : public vcg::tri::TriEdgeCollapseQuadric<
  DecMesh, DecVertexPair, DecCollapse, vcg::tri::QInfoStandard<DecVertex>
> {
public:
  typedef vcg::tri::TriEdgeCollapseQuadric<
    DecMesh, DecVertexPair, DecCollapse, vcg::tri::QInfoStandard<DecVertex>
  > TECQ;
  inline DecCollapse(const DecVertexPair &p, int i, vcg::BaseParameterClass *pp)
    : TECQ(p, i, pp) {}
};

} // namespace

// [[Rcpp::export]]
Rcpp::List vcgDecimate(
    SEXP vb_, SEXP it_, int targetFaces, bool preserveTopology,
    bool preserveBoundary, bool normalCheck, double qualityThreshold,
    bool verbose)
{
  try {
    DecMesh m;
    int check = ravetools::IOMesh<DecMesh>::vcgReadR(m, vb_, it_);
    if (check != 0) {
      Rcpp::stop("vcg_decimate: mesh has no faces and/or no vertices");
    }

    vcg::tri::Clean<DecMesh>::RemoveDuplicateVertex(m);
    vcg::tri::Clean<DecMesh>::RemoveUnreferencedVertex(m);
    vcg::tri::Allocator<DecMesh>::CompactEveryVector(m);
    vcg::tri::UpdateBounding<DecMesh>::Box(m);

    vcg::tri::TriEdgeCollapseQuadricParameter params;
    params.PreserveTopology = preserveTopology;
    params.PreserveBoundary = preserveBoundary;
    params.NormalCheck = normalCheck;
    params.QualityCheck = true;
    params.QualityThr = qualityThreshold;
    params.OptimalPlacement = true;
    params.ScaleIndependent = true;

    if (verbose) {
      Rprintf("vcg_decimate: %d vertices, %d faces; target %d faces\n",
              m.vn, m.fn, targetFaces);
    }

    vcg::LocalOptimization<DecMesh> session(m, &params);
    session.Init<DecCollapse>();
    // stale candidates stay on the heap until it outgrows HeapSimplexRatio
    // times the face count; lowering vcglib's 4 to 2 purges them twice as
    // often, which cut the peak memory of a 3.2 M vertex mesh by ~0.5 GB
    // and did not slow the decimation
    session.HeapSimplexRatio = 2.0f;
    session.SetTargetSimplices(targetFaces);
    session.SetTimeBudget(0.5f);
    while (m.fn > targetFaces) {
      const int before = m.fn;
      session.DoOptimization();
      Rcpp::checkUserInterrupt();
      // the heap can run dry before the target when the constraints (topology,
      // boundary, normals) forbid every remaining collapse
      if (m.fn >= before) break;
    }
    // restores the write flags of preserved boundary vertices; the static list
    // vcglib keeps them in would otherwise point into this mesh after it is gone
    session.Finalize<DecCollapse>();
    DecCollapse::WV().clear();

    vcg::tri::Allocator<DecMesh>::CompactVertexVector(m);
    vcg::tri::Allocator<DecMesh>::CompactFaceVector(m);
    vcg::tri::UpdateNormal<DecMesh>::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal<DecMesh>::NormalizePerVertex(m);

    if (verbose) {
      Rprintf("vcg_decimate: result %d vertices, %d faces (estimated error %g)\n",
              m.vn, m.fn, session.currMetric);
    }

    return ravetools::IOMesh<DecMesh>::vcgToR(m, true);

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return Rcpp::List(); // -Wall
}
