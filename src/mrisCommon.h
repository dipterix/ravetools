#ifndef RAVETOOLS_MRIS_COMMON_H
#define RAVETOOLS_MRIS_COMMON_H

// mrisCommon.h
//
// Shared mesh representation and helpers for the cortical-surface algorithm
// implementations (mris_inflate, mris_smooth, mris_sphere, ...).
//
// These are *not* built on VCG: the underlying surface-relaxation algorithms
// operate on a flat structure-of-arrays representation with 1-/2-ring vertex
// adjacency, which is simpler and faster to build directly from the face list
// than to route through VCG's general-purpose mesh machinery. Sharing this
// module across the mris_* ports avoids re-deriving the same adjacency,
// normal, area, and integration bookkeeping in each one.
//
// Reference: Fischl B, Sereno MI, Dale AM (1999). Cortical surface-based
// analysis II: Inflation, flattening, and a surface-based coordinate system.
// NeuroImage, 9(2), 195-207.

#include <Rcpp.h>
#include <vector>
#include <cmath>

namespace ravetools {

// ---------------------------------------------------------------------------
// Flat SOA surface representation (the subset of per-vertex/per-face state
// needed by the mris_* ports).
// ---------------------------------------------------------------------------
struct MrisMesh {
    int nv = 0, nf = 0;

    // Current vertex positions
    std::vector<float> x, y, z;

    // Face indices (0-based)
    std::vector<int> f0, f1, f2;

    // Per-vertex normals (kept current via mrisComputeNormals)
    std::vector<float> nx, ny, nz;

    // Gradient accumulator (cleared at the start of each iteration)
    std::vector<float> dx, dy, dz;

    // Momentum velocity accumulator for the momentum-integration time step
    std::vector<float> odx, ody, odz;

    // 1-ring adjacency: direct edge-sharing neighbors.
    std::vector<std::vector<int>> nbrs1;

    // 2-ring adjacency: vertices within 2 hops, used (along with avg_nbrs)
    // to normalize and smooth the per-vertex gradient terms.
    std::vector<std::vector<int>> nbrs2;

    // Distances to 2-ring neighbors: "orig" is the reference the distance
    // term tries to preserve, "curr" is recomputed every iteration.
    std::vector<std::vector<float>> dist_orig;
    std::vector<std::vector<float>> dist_curr;

    // Area statistics (recomputed via mrisComputeArea)
    float orig_area  = 0.0f;   // reference area, set once, never changes
    float total_area = 0.0f;
    float neg_area   = 0.0f;

    // Average 2-ring neighborhood size; used to normalize the distance-term
    // gradient.
    float avg_nbrs = 0.0f;
};

// ---------------------------------------------------------------------------
// Construction / validation
// ---------------------------------------------------------------------------

// Load vertex positions and faces from R matrices (vb_: 3 x nv numeric,
// it_: 3 x nf integer, already 0-based, caller subtracts 1 in R). Allocates
// and zeroes all derived state vectors (normals, gradient, momentum).
void mrisLoadMesh(MrisMesh &m, SEXP vb_, SEXP it_);

// Validate that the mesh is closed and manifold (every edge shared by exactly
// two faces). `caller` names the function in the error message (e.g.
// "mris_inflate"). See the detailed rationale in mrisCommon.cpp: the
// underlying surface-relaxation algorithms are only intended to run on
// topologically-correct, closed (watertight), genus-0 surfaces and are not
// robust to boundaries, handles, or non-manifold edges.
void mrisCheckClosedManifold(const MrisMesh &m, const char *caller);

// ---------------------------------------------------------------------------
// Per-face signed area and outward unit normal
// ---------------------------------------------------------------------------

struct FaceGeometry { float nx, ny, nz, area; };

// Per-face signed area and outward-pointing unit normal, oriented for closed
// ellipsoid/sphere-status surfaces (the only kind the mris_* ports operate
// on): a face's raw cross-product normal n = (v1-v0) x (v2-v0) is treated as
// pointing inward, and its area negated, its normal flipped, when
// dot(v0+v1+v2, n) < 0, i.e. when n points towards rather than away from the
// surface's center. This requires the mesh to be centered at the origin
// first (mrisCenter).
void mrisComputeFaceGeometry(const MrisMesh &m, std::vector<FaceGeometry> &fg);

// ---------------------------------------------------------------------------
// Topology / metric properties
// ---------------------------------------------------------------------------

// Build 1-ring adjacency (nbrs1). If `two_ring` is true, also builds the
// 2-ring adjacency (nbrs2) and computes avg_nbrs.
void mrisBuildAdjacency(MrisMesh &m, bool two_ring);

// Per-vertex normals: area-weighted average of surrounding face normals,
// normalized.
void mrisComputeNormals(MrisMesh &m);

// Total and "negative" mesh area, via mrisComputeFaceGeometry's signed
// per-face areas: total_area sums the positive-area faces, neg_area sums the
// magnitude of the negative-area (folded/inverted) ones; the reported total
// surface area excludes the negative-area faces.
void mrisComputeArea(MrisMesh &m);

// Current Euclidean distances from each vertex to its 2-ring neighbors
// (requires nbrs2; fills dist_curr).
void mrisComputeCurrentDistances(MrisMesh &m);

// Populate dist_orig from the current vertex positions (snapshot taken once,
// after initial scaling).
void mrisStoreOriginalDistances(MrisMesh &m);

// ---------------------------------------------------------------------------
// Geometric transforms
// ---------------------------------------------------------------------------

// Uniform 3-D scale about the origin so that total_area == target_area.
void mrisScaleToArea(MrisMesh &m, float target_area);

// Translate the centroid to the origin.
void mrisCenter(MrisMesh &m);

// Radially project every vertex onto a sphere of radius `r` centered at the
// origin: x <- x * r / |x| (the caller is responsible for centering first).
void mrisProjectOntoSphere(MrisMesh &m, float r);

// ---------------------------------------------------------------------------
// Iterative-integration building blocks (shared by mris_inflate / mris_sphere)
// ---------------------------------------------------------------------------

// In-place include-self 1-ring averaging of a vector field, `n_averages`
// passes:  new[v] = (field[v] + sum(field[n] for n in ring1[v])) / (1 + |ring1[v]|)
// This is the shared smoothing/averaging step used both for the force
// gradient (mris_inflate / mris_sphere) and for vertex positions
// (mris_smooth).
void mrisAverageField(const MrisMesh &m, std::vector<float> &fx,
                      std::vector<float> &fy, std::vector<float> &fz,
                      int n_averages);

// Apply momentum to the gradient accumulator (dx/dy/dz), clamp the per-vertex
// displacement to `max_disp_mm`, and move the vertices (both mris_inflate
// and the reduced mris_sphere unfolding use this as the final step of each
// iteration, with no further gradient averaging applied here).
void mrisMomentumStep(MrisMesh &m, float momentum, float dt, float max_disp_mm);

// RMS of each vertex's height above the tangent plane defined by its 1-ring
// neighbor centroid and normal; used as an early-stopping criterion.
float mrisRmsTPHeight(const MrisMesh &m);

// ---------------------------------------------------------------------------
// Output packing
// ---------------------------------------------------------------------------

// Pack vertex positions and normals into 3 x nv numeric matrices, returned to
// R as `vb` and `normals` (the R wrapper re-attaches the 1-based face matrix
// and the homogeneous coordinate row).
Rcpp::List mrisPackMeshOutput(const MrisMesh &m);

} // namespace ravetools

#endif // RAVETOOLS_MRIS_COMMON_H
