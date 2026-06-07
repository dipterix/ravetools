// mris.cpp
//
// Cortical-surface algorithm ports built on the shared representation and
// helpers in mrisCommon.h: surface smoothing (mris_smooth), cortical surface
// inflation (mris_inflate), and spherical projection / metric-distortion
// relaxation (mris_sphere). Grouped into a single translation unit,the
// same way the VCG-based routines are grouped in vcgCommon.cpp,so the
// shared header and integration machinery are only compiled once, which
// keeps the build smaller and faster.
//
// Reference: Fischl B, Sereno MI, Dale AM (1999). Cortical surface-based
// analysis II: Inflation, flattening, and a surface-based coordinate system.
// NeuroImage, 9(2), 195-207.

#include "mrisCommon.h"
#include <algorithm>

using namespace ravetools;

namespace {

// Per-vertex momentum clamp shared by mris_inflate and mris_sphere (limits
// how far a vertex can move in a single integration step).
const float MAX_MOMENTUM_MM = 1.0f;

// ===========================================================================
// mris_inflate helpers
// ===========================================================================

// Canonical brain surface area used as the inflation target (110 cm^2)
const float DEFAULT_BRAIN_AREA = 110000.0f;

// ---------------------------------------------------------------------------
// FORCE TERM 1: distance-preserving restoring force.
//
// For each vertex v, for all 2-ring neighbors n:
//   d0    = dist_orig[n] / scale    (desired distance, area-corrected)
//   dt    = dist_curr[n]            (actual distance)
//   delta = dt - d0
//   grad += delta * normalize(n - v)
// grad /= avg_nbrs                  (normalize by average neighborhood size)
// Remove normal component (tangential only).
// dx += l_dist * grad
//
// l_dist is passed pre-scaled by sqrt(n_averages) by the caller.
// ---------------------------------------------------------------------------
void inflateDistanceTerm(MrisMesh &m, float l_dist)
{
    if (l_dist == 0.0f) return;

    // area-correction scale: discount folded/negative-area faces when present
    float scale = (m.neg_area < m.total_area && m.total_area > 0.0f)
                    ? std::sqrt(m.orig_area / (m.total_area - m.neg_area))
                    : std::sqrt(m.orig_area / m.total_area);

    float norm = (m.avg_nbrs > 0.0f) ? 1.0f / m.avg_nbrs : 1.0f;

    for (int vi = 0; vi < m.nv; vi++) {
        const auto &nb    = m.nbrs2[vi];
        const auto &d_orig = m.dist_orig[vi];
        const auto &d_curr = m.dist_curr[vi];
        int nn = (int)nb.size();

        float fdx = 0.0f, fdy = 0.0f, fdz = 0.0f;
        for (int ni = 0; ni < nn; ni++) {
            int nj = nb[ni];
            float d0    = d_orig[ni] / scale;
            float delta = d_curr[ni] - d0;

            float yx = m.x[nj]-m.x[vi], yy = m.y[nj]-m.y[vi], yz = m.z[nj]-m.z[vi];
            float yl = std::sqrt(yx*yx + yy*yy + yz*yz);
            if (yl < 1e-10f) continue;
            yx /= yl; yy /= yl; yz /= yl;

            fdx += delta * yx;
            fdy += delta * yy;
            fdz += delta * yz;
        }
        fdx *= norm; fdy *= norm; fdz *= norm;

        // Remove normal component (keep only tangential force)
        float nc = fdx*m.nx[vi] + fdy*m.ny[vi] + fdz*m.nz[vi];
        fdx -= nc*m.nx[vi]; fdy -= nc*m.ny[vi]; fdz -= nc*m.nz[vi];

        m.dx[vi] += l_dist * fdx;
        m.dy[vi] += l_dist * fdy;
        m.dz[vi] += l_dist * fdz;
    }
}

// ---------------------------------------------------------------------------
// FORCE TERM 2: normalized Laplacian spring term.
//
// Added AFTER mrisAverageField(dx,dy,dz) is applied to the distance gradient,
// so this term is not smoothed (acts locally).
// For each vertex v with 1-ring neighbors:
//   spring_disp = dist_scale/N * sum(n - v)   (scaled Laplacian displacement)
//   dx += l_spring_norm * spring_disp
// where dist_scale = sqrt(orig_area / total_area).
//
// CRITICAL second pass (easy to miss): the mesh-wide average of the
// normal-component of this term, dot_avg, is computed and then SUBTRACTED
// from every vertex's gradient (dx -= dot_avg * normal). This removes the net
// radial bias that a Laplacian-type spring term otherwise introduces (it both
// smooths AND uniformly contracts/expands the surface along normals). Without
// this correction the bias accumulates under momentum integration and
// produces runaway, direction-dependent distortion.
// ---------------------------------------------------------------------------
void normalizedSpringTerm(MrisMesh &m, float l_spring_norm)
{
    if (l_spring_norm == 0.0f) return;
    float dist_scale = std::sqrt(m.orig_area / m.total_area);

    std::vector<float> sxv(m.nv, 0.0f), syv(m.nv, 0.0f), szv(m.nv, 0.0f);
    double dot_total = 0.0;

    for (int vi = 0; vi < m.nv; vi++) {
        const auto &nb = m.nbrs1[vi];
        int nn = (int)nb.size();
        if (nn == 0) continue;
        float sx = 0.0f, sy = 0.0f, sz = 0.0f;
        for (int nj : nb) { sx += m.x[nj]-m.x[vi]; sy += m.y[nj]-m.y[vi]; sz += m.z[nj]-m.z[vi]; }
        float mul = dist_scale / nn;
        sx *= mul; sy *= mul; sz *= mul;
        sxv[vi] = sx; syv[vi] = sy; szv[vi] = sz;

        dot_total += (double)(l_spring_norm * (m.nx[vi]*sx + m.ny[vi]*sy + m.nz[vi]*sz));

        m.dx[vi] += l_spring_norm * sx;
        m.dy[vi] += l_spring_norm * sy;
        m.dz[vi] += l_spring_norm * sz;
    }

    float dot_avg = (m.nv > 0) ? (float)(dot_total / m.nv) : 0.0f;

    for (int vi = 0; vi < m.nv; vi++) {
        m.dx[vi] -= dot_avg * m.nx[vi];
        m.dy[vi] -= dot_avg * m.ny[vi];
        m.dz[vi] -= dot_avg * m.nz[vi];
    }
}

// ---------------------------------------------------------------------------
// SULCAL DEPTH
//
// Accumulates the normal-projected step velocity into the sulc field.
// Uses normals from BEFORE the current position update (stale at call time):
// move vertices -> accumulate depth -> refresh metric properties.
// ---------------------------------------------------------------------------
void trackSulcalDepth(MrisMesh &m, std::vector<float> &sulc)
{
    for (int vi = 0; vi < m.nv; vi++) {
        sulc[vi] += m.odx[vi]*m.nx[vi] + m.ody[vi]*m.ny[vi] + m.odz[vi]*m.nz[vi];
    }
}

// ===========================================================================
// mris_sphere helpers
// ===========================================================================

// Canonical sphere radius, the default target_radius in mris_sphere.
const float DEFAULT_RADIUS = 100.0f;

// ---------------------------------------------------------------------------
// FORCE TERM 1: distance-preserving restoring force (sphere-status branch).
//
// Identical in structure to the distance term in mris_inflate, EXCEPT for the
// area-correction `scale`: the sphere-status branch always uses the plain
// sqrt(orig_area / total_area),it does not take the neg_area-aware branch
// that the generic (non-sphere) path uses.
// ---------------------------------------------------------------------------
void sphereDistanceTerm(MrisMesh &m, float l_dist)
{
    if (l_dist == 0.0f) return;

    float scale = (m.total_area > 0.0f) ? std::sqrt(m.orig_area / m.total_area) : 1.0f;
    float norm  = (m.avg_nbrs > 0.0f) ? 1.0f / m.avg_nbrs : 1.0f;

    for (int vi = 0; vi < m.nv; vi++) {
        const auto &nb     = m.nbrs2[vi];
        const auto &d_orig = m.dist_orig[vi];
        const auto &d_curr = m.dist_curr[vi];
        int nn = (int)nb.size();

        float fdx = 0.0f, fdy = 0.0f, fdz = 0.0f;
        for (int ni = 0; ni < nn; ni++) {
            int nj = nb[ni];
            float d0    = d_orig[ni] / scale;
            float delta = d_curr[ni] - d0;

            float yx = m.x[nj]-m.x[vi], yy = m.y[nj]-m.y[vi], yz = m.z[nj]-m.z[vi];
            float yl = std::sqrt(yx*yx + yy*yy + yz*yz);
            if (yl < 1e-10f) continue;
            yx /= yl; yy /= yl; yz /= yl;

            fdx += delta * yx;
            fdy += delta * yy;
            fdz += delta * yz;
        }
        fdx *= norm; fdy *= norm; fdz *= norm;

        // Remove normal component (keep only tangential force)
        float nc = fdx*m.nx[vi] + fdy*m.ny[vi] + fdz*m.nz[vi];
        fdx -= nc*m.nx[vi]; fdy -= nc*m.ny[vi]; fdz -= nc*m.nz[vi];

        m.dx[vi] += l_dist * fdx;
        m.dy[vi] += l_dist * fdy;
        m.dz[vi] += l_dist * fdz;
    }
}

// ---------------------------------------------------------------------------
// FORCE TERM 2: repulsive force on folded/negative-area faces only (the
// l_area-only branch of the energy functional; l_angle and l_parea are zero
// for this reduced unfolding pass).
//
// This term is a *repulsive* force that acts only on faces whose (rescaled)
// area has gone non-positive (i.e. folded/inverted faces): it pushes their
// vertices apart proportionally to how negative the face has become, which is
// exactly the "unfolding" mechanism. Faces with positive area contribute zero
// (delta stays 0,see the `if (area <= 0.0f)` guard below).
//
// For each face f = (v0, v1, v2) with edges a = v1-v0, b = v2-v0 and unit
// normal n:
//   area_scale = orig_area_total / total_area      (METRIC_SCALE, sphere status)
//   area       = area_scale * face_area
//   delta      = l_area * (area - orig_face_area)  [only if area <= 0]
//   v0 += delta * (-(a x n) + (b x n))
//   v1 += delta * (-(b x n))
//   v2 += delta *   (a x n)
// ---------------------------------------------------------------------------
void areaTerm(MrisMesh &m, const std::vector<FaceGeometry> &orig_fg,
              const std::vector<FaceGeometry> &fg, float l_area)
{
    if (l_area == 0.0f) return;

    float area_scale = (m.total_area > 0.0f) ? (m.orig_area / m.total_area) : 1.0f;

    for (int fi = 0; fi < m.nf; fi++) {
        float area = area_scale * fg[fi].area;
        if (area > 0.0f) continue;   // only folded/negative faces are pushed

        float delta = l_area * (area - orig_fg[fi].area);

        int i0 = m.f0[fi], i1 = m.f1[fi], i2 = m.f2[fi];
        float ax = m.x[i1]-m.x[i0], ay = m.y[i1]-m.y[i0], az = m.z[i1]-m.z[i0];
        float bx = m.x[i2]-m.x[i0], by = m.y[i2]-m.y[i0], bz = m.z[i2]-m.z[i0];
        float nx = fg[fi].nx, ny = fg[fi].ny, nz = fg[fi].nz;

        float axnx = ay*nz - az*ny, axny = az*nx - ax*nz, axnz = ax*ny - ay*nx;
        float bxnx = by*nz - bz*ny, bxny = bz*nx - bx*nz, bxnz = bx*ny - by*nx;

        m.dx[i0] += delta * (-axnx + bxnx);
        m.dy[i0] += delta * (-axny + bxny);
        m.dz[i0] += delta * (-axnz + bxnz);

        m.dx[i1] += delta * (-bxnx);
        m.dy[i1] += delta * (-bxny);
        m.dz[i1] += delta * (-bxnz);

        m.dx[i2] += delta * axnx;
        m.dy[i2] += delta * axny;
        m.dz[i2] += delta * axnz;
    }
}

} // namespace

// ===========================================================================
// mris_smooth
//
// Surface tessellation smoothing. For each of `npasses` outer passes:
//
//   average vertex positions over `niterations` rounds of include-self
//   1-ring averaging (the per-vertex mean of itself and its directly-
//   connected neighbors);
//   recompute normals and surface-area statistics;
//   if `rescale`, uniformly scale back to the original surface area.
// ===========================================================================

// ---------------------------------------------------------------------------
// Matches the default mris_smooth call: niterations=10, npasses=1, rescale=FALSE
//
// Returns a list with:
//   vb      - 3 x nv smoothed vertex positions
//   normals - 3 x nv vertex normals of the smoothed surface
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mrisSmooth(
    SEXP        vb_,
    SEXP        it_,
    int         niterations = 10,
    int         npasses     = 1,
    bool        rescale     = false,
    bool        verbose     = false
) {
    MrisMesh m;
    mrisLoadMesh(m, vb_, it_);

    mrisCheckClosedManifold(m, "mris_smooth");

    if (verbose) Rprintf("mris_smooth: building adjacency for %d vertices, %d faces\n", m.nv, m.nf);
    mrisBuildAdjacency(m, /* two_ring = */ false);

    mrisComputeNormals(m);
    mrisComputeArea(m);
    m.orig_area = m.total_area;   // reference area, used to rescale back if requested

    for (int pass = 0; pass < npasses; pass++) {
        Rcpp::checkUserInterrupt();

        if (verbose) Rprintf("mris_smooth: pass %d/%d: averaging vertex positions over %d iterations\n",
                             pass + 1, npasses, niterations);

        // `niterations` rounds of include-self 1-ring averaging applied
        // directly to the vertex coordinates (the same operation mris_inflate
        // applies to its force gradient, generalized here as mrisAverageField).
        mrisAverageField(m, m.x, m.y, m.z, niterations);

        // Refresh normals and surface-area statistics for the smoothed mesh
        mrisComputeNormals(m);
        mrisComputeArea(m);

        // Restore the original surface area
        if (rescale) {
            mrisScaleToArea(m, m.orig_area);
            mrisComputeNormals(m);
        }
    }

    if (verbose) Rprintf("mris_smooth: done\n");

    return mrisPackMeshOutput(m);
}

// ===========================================================================
// mris_inflate
//
// Cortical surface inflation: iteratively relaxes a closed cortical-surface
// mesh into a smoother, more compact shape while tracking how far each
// vertex moves inward along the surface normal.
// ===========================================================================

// ---------------------------------------------------------------------------
// Main exported function
//
// Matches the default mris_inflate call:
//   n_averages=16, niterations=10, l_spring_norm=1.0, l_dist=0.1,
//   momentum=0.9, dt=0.9, desired_rms=0.015, scale_brain=true
//
// Returns a list with:
//   vb      - 3 x nv inflated vertex positions
//   normals - 3 x nv vertex normals of inflated surface
//   sulc    - nv sulcal depth values (zero-mean, in mesh units)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mrisInflate(
    SEXP        vb_,
    SEXP        it_,
    int         n_averages    = 16,
    int         niterations   = 10,
    double      l_spring_norm = 1.0,
    double      l_dist        = 0.1,
    double      momentum      = 0.9,
    double      dt            = 0.9,
    double      desired_rms   = 0.015,
    bool        scale_brain   = true,
    bool        verbose       = false
) {
    MrisMesh m;
    mrisLoadMesh(m, vb_, it_);

    // --- validate topology (closed, manifold, genus-0,see comment in
    //     mrisCheckClosedManifold for why this is a hard precondition) ---
    mrisCheckClosedManifold(m, "mris_inflate");

    std::vector<float> sulc(m.nv, 0.0f);   // per-vertex sulcal depth accumulator

    // --- Pre-processing ---

    if (verbose) Rprintf("mris_inflate: building adjacency for %d vertices, %d faces\n", m.nv, m.nf);
    mrisBuildAdjacency(m, /* two_ring = */ true);

    // Normals and area from the input mesh
    mrisComputeNormals(m);
    mrisComputeArea(m);

    // Uniform scale so total_area == target_area.
    // We store orig_area as the post-scale area so that dist_scale starts at 1.
    float target_area = scale_brain ? DEFAULT_BRAIN_AREA : m.total_area;
    mrisScaleToArea(m, target_area);
    m.orig_area = target_area;

    // Normals don't change direction under uniform scaling, but recompute
    // to ensure they are consistent with the updated positions.
    mrisComputeNormals(m);

    // Reference distances from the scaled positions (these are what l_dist
    // tries to preserve).
    if (verbose) Rprintf("mris_inflate: computing original distances\n");
    mrisStoreOriginalDistances(m);

    // --- Inflation loop ---
    //
    // Outer loop: n_avg from n_averages down to 0, halved each level.
    // Inner loop: niterations per level (may break early on RMS criterion).
    //
    // l_dist is scaled by sqrt(n_avg) per level. At n_avg=0 only the spring
    // term acts.

    float l_dist_base        = (float)l_dist;
    float l_spring_norm_f    = (float)l_spring_norm;
    float momentum_f         = (float)momentum;
    float dt_f               = (float)dt;
    float desired_rms_f      = (float)desired_rms;
    int   total_iters        = 0;

    if (verbose) {
        float rms0 = mrisRmsTPHeight(m);
        Rprintf("mris_inflate: starting inflation, initial RMS=%.4f (target=%.4f)\n",
                rms0, desired_rms_f);
    }

    int n_avg = n_averages;
    while (true) {
        float l_dist_eff = l_dist_base * std::sqrt((float)n_avg);

        for (int n = 0; n < niterations; n++) {
            Rcpp::checkUserInterrupt();

            // 1. Clear gradient
            std::fill(m.dx.begin(), m.dx.end(), 0.0f);
            std::fill(m.dy.begin(), m.dy.end(), 0.0f);
            std::fill(m.dz.begin(), m.dz.end(), 0.0f);

            // 2. Distance term (metric-preserving restoring force)
            inflateDistanceTerm(m, l_dist_eff);

            // 3. Smooth the distance gradient over n_avg passes of 1-ring averaging.
            //    The spring term below is added AFTER this and is NOT averaged.
            mrisAverageField(m, m.dx, m.dy, m.dz, n_avg);

            // 4. Normalized spring term (local smoothing, not averaged)
            normalizedSpringTerm(m, l_spring_norm_f);

            // 5. Apply momentum and move vertices
            mrisMomentumStep(m, momentum_f, dt_f, MAX_MOMENTUM_MM);

            // 6. Accumulate sulcal depth (uses pre-update normals, before recompute)
            trackSulcalDepth(m, sulc);

            // 7. Update metric properties for next iteration
            mrisComputeCurrentDistances(m);
            mrisComputeArea(m);
            mrisComputeNormals(m);

            total_iters++;

            float rms = mrisRmsTPHeight(m);
            if (verbose && (total_iters % 5 == 0)) {
                Rprintf("  iter %3d (navg=%2d, ldist=%.3f): RMS=%.5f\n",
                        total_iters, n_avg, l_dist_eff, rms);
            }
            if (rms < desired_rms_f) break;   // inner-loop early-stopping criterion only
        }

        if (n_avg == 0) break;
        n_avg /= 2;
    }

    if (verbose) Rprintf("mris_inflate: done (%d total iterations)\n", total_iters);

    // --- Post-processing ---

    // Translate centroid to origin
    mrisCenter(m);

    // Restore the canonical brain surface area
    mrisComputeArea(m);
    mrisScaleToArea(m, m.orig_area);

    // Recompute normals for the returned surface
    mrisComputeNormals(m);

    // Zero-mean the sulcal depth map
    double sulc_mean = 0.0;
    for (int vi = 0; vi < m.nv; vi++) sulc_mean += sulc[vi];
    sulc_mean /= m.nv;
    for (int vi = 0; vi < m.nv; vi++) sulc[vi] -= (float)sulc_mean;

    // --- Pack output ---
    Rcpp::List out = mrisPackMeshOutput(m);
    Rcpp::NumericVector out_sulc(m.nv);
    for (int vi = 0; vi < m.nv; vi++) out_sulc[vi] = sulc[vi];
    out["sulc"] = out_sulc;

    return out;
}

// ===========================================================================
// mris_sphere
//
// Spherical projection and metric-distortion relaxation ("unfolding") of a
// surface mesh.
//
// This is a *reduced* implementation, not a faithful reproduction of any
// particular reference pipeline. The radial projection step is implemented
// faithfully (see mrisProjectOntoSphere / sphereDistanceTerm above); the full
// unfolding procedure described in the literature drives a multi-resolution,
// line-minimization optimization through several thousand lines of machinery
// (seven weighted energy terms, l_dist, l_area, l_angle, l_neg, l_spring,
// l_curv, l_nlarea, against a separate "original"/white-matter reference
// surface, with multi-resolution remeshing, ...). Reproducing that faithfully
// is out of scope for this package. Instead this implements an unfolding loop
// that keeps only the two dominant default energy terms:
//   l_dist = 1.0  (sphereDistanceTerm - metric-preserving restoring force)
//   l_area = 1.0  (areaTerm           - repulsive force on folded/negative-
//                                         area faces only)
// and integrates them with the SAME momentum-integration machinery as
// mris_inflate (gradient averaging, momentum, clamped displacement), the
// same machinery the reference pipeline also uses for the inflation phase
// that normally precedes unfolding.
//
// A further simplification: the reference pipeline's distance term targets a
// separately-loaded *white-matter* reference surface; since this package's
// input is a single mesh, the input mesh's own metric is used as the
// reference instead, exactly the convention mris_inflate already uses (and
// the practical analogue, since the white-matter surface is normally what
// gets inflated to produce the input to this kind of unfolding step).
// ===========================================================================

// ---------------------------------------------------------------------------
// Main exported function.
//
// Pipeline (see section header above for the rationale behind each
// simplification):
//   1. Center, snapshot the input mesh's own metric as the distance-term
//      reference (dist_orig) and per-face areas as the area-term reference.
//   2. Project radially onto a sphere of radius `target_radius` (faithful
//      radial projection, see mrisProjectOntoSphere).
//   3. Unfolding loop, `niterations` rounds of:
//        clear gradient -> distance term -> average gradient (n_averages
//        passes of 1-ring smoothing, shared with mris_inflate's gradient
//        averaging step) -> area term (added after averaging, acts locally,
//        mirroring how mris_inflate adds its spring term) -> momentum step ->
//        refresh metric properties.
//
// Defaults mirror the reference pipeline's *unfolding-pass* parameters
// (l_dist=1.0, l_area=1.0, niterations=25, momentum=0.9, dt=0.05,
// target_radius=100), with n_averages reduced from the production default of
// 1024 (tuned for ~150k-vertex cortical surfaces) to a size more practical
// for typical 'mesh3d' inputs, documented here as a deliberate departure,
// not an oversight.
//
// Returns a list with:
//   vb      - 3 x nv spherical vertex positions
//   normals - 3 x nv vertex normals of the resulting sphere
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mrisSphere(
    SEXP        vb_,
    SEXP        it_,
    double      target_radius = 100.0,
    int         n_averages    = 64,
    int         niterations   = 25,
    double      l_dist        = 1.0,
    double      l_area        = 1.0,
    double      momentum      = 0.9,
    double      dt            = 0.05,
    bool        verbose       = false
) {
    MrisMesh m;
    mrisLoadMesh(m, vb_, it_);

    mrisCheckClosedManifold(m, "mris_sphere");

    if (verbose) Rprintf("mris_sphere: building adjacency for %d vertices, %d faces\n", m.nv, m.nf);
    mrisBuildAdjacency(m, /* two_ring = */ true);

    mrisComputeNormals(m);
    mrisComputeArea(m);
    mrisCenter(m);

    // --- Reference metric: the input mesh's own geometry (see section header
    //     above for why this stands in for a separately-loaded white-matter
    //     reference surface). ---
    m.orig_area = m.total_area;
    mrisStoreOriginalDistances(m);

    std::vector<FaceGeometry> orig_fg;
    mrisComputeFaceGeometry(m, orig_fg);

    float radius = (float)target_radius;
    if (radius <= 0.0f) radius = DEFAULT_RADIUS;

    // --- Faithful radial projection onto the target sphere ---
    if (verbose) Rprintf("mris_sphere: projecting onto sphere of radius %.2f\n", radius);
    mrisProjectOntoSphere(m, radius);

    mrisComputeNormals(m);
    mrisComputeArea(m);
    mrisComputeCurrentDistances(m);

    // --- Unfolding loop ---

    float l_dist_f   = (float)l_dist;
    float l_area_f   = (float)l_area;
    float momentum_f = (float)momentum;
    float dt_f       = (float)dt;

    if (verbose) {
        Rprintf("mris_sphere: starting unfolding, initial neg_area=%.4f (total_area=%.4f)\n",
                m.neg_area, m.total_area);
    }

    std::vector<FaceGeometry> fg;
    for (int iter = 0; iter < niterations; iter++) {
        Rcpp::checkUserInterrupt();

        // 1. Clear gradient
        std::fill(m.dx.begin(), m.dx.end(), 0.0f);
        std::fill(m.dy.begin(), m.dy.end(), 0.0f);
        std::fill(m.dz.begin(), m.dz.end(), 0.0f);

        // 2. Distance term (metric-preserving restoring force)
        sphereDistanceTerm(m, l_dist_f);

        // 3. Smooth the distance gradient (mirrors mris_inflate: the area
        //    term below is added AFTER averaging and is NOT smoothed, just
        //    like mris_inflate's spring term).
        mrisAverageField(m, m.dx, m.dy, m.dz, n_averages);

        // 4. Area term (local repulsion on folded faces, not averaged)
        mrisComputeFaceGeometry(m, fg);
        areaTerm(m, orig_fg, fg, l_area_f);

        // 5. Apply momentum and move vertices
        mrisMomentumStep(m, momentum_f, dt_f, MAX_MOMENTUM_MM);

        // 6. Update metric properties for next iteration
        mrisComputeCurrentDistances(m);
        mrisComputeArea(m);
        mrisComputeNormals(m);

        if (verbose && ((iter + 1) % 5 == 0)) {
            Rprintf("  iter %3d: neg_area=%.5f (total_area=%.4f)\n",
                    iter + 1, m.neg_area, m.total_area);
        }
    }

    if (verbose) Rprintf("mris_sphere: done (%d iterations)\n", niterations);

    return mrisPackMeshOutput(m);
}
