// mris.cpp
//
// Cortical-surface algorithm ports built on the shared representation and
// helpers in mrisCommon.h: surface smoothing (mris_smooth), cortical surface
// inflation (mris_inflate), and spherical projection / metric-distortion
// relaxation (mris_sphere). Grouped into a single translation unit, the same
// way the VCG-based routines are grouped in vcgCommon.cpp, so the shared
// header and integration machinery are only compiled once, which keeps the
// build smaller and faster.
//
// Reference: Fischl B, Sereno MI, Dale AM (1999). Cortical surface-based
// analysis II: Inflation, flattening, and a surface-based coordinate system.
// NeuroImage, 9(2), 195-207.

#include "mrisCommon.h"
#include <algorithm>
#include <unordered_map>
#include <cstdint>

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
// sqrt(orig_area / total_area); it does not take the neg_area-aware branch
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

// ===========================================================================
// mris_make_surfaces helpers
// ===========================================================================

// Per-vertex displacement clamp for the intensity-driven placement passes;
// reuses the same 1 mm per-step cap as mris_inflate / mris_sphere.
const float MAX_SURFACE_STEP_MM = MAX_MOMENTUM_MM;

// ---------------------------------------------------------------------------
// Trilinear sampler for a 3-D intensity volume in surface ('RAS') space.
//
// `r2i` is the affine that maps a surface-space (x, y, z) point to its
// fractional volume-index ('IJK', 0-based, column-major) coordinates; it is
// the inverse of the volume's IJK-to-RAS transform, flattened row-major (only
// the first three rows are needed, the fourth is always [0 0 0 1]):
//   i = r2i[0]*x + r2i[1]*y + r2i[2]*z  + r2i[3]
//   j = r2i[4]*x + r2i[5]*y + r2i[6]*z  + r2i[7]
//   k = r2i[8]*x + r2i[9]*y + r2i[10]*z + r2i[11]
// ---------------------------------------------------------------------------
struct IntensityVolume {
    const double *data = nullptr;
    int nx = 0, ny = 0, nz = 0;
    float r2i[12] = {0.0f};

    inline bool sample(float x, float y, float z, float &out) const
    {
        float fi = r2i[0]*x + r2i[1]*y + r2i[2]*z  + r2i[3];
        float fj = r2i[4]*x + r2i[5]*y + r2i[6]*z  + r2i[7];
        float fk = r2i[8]*x + r2i[9]*y + r2i[10]*z + r2i[11];

        int i0 = (int)std::floor(fi), j0 = (int)std::floor(fj), k0 = (int)std::floor(fk);
        if (i0 < 0 || j0 < 0 || k0 < 0 || i0 + 1 >= nx || j0 + 1 >= ny || k0 + 1 >= nz) return false;

        float di = fi - i0, dj = fj - j0, dk = fk - k0;

        auto at = [&](int i, int j, int k) -> double {
            return data[(size_t)i + (size_t)nx * ((size_t)j + (size_t)ny * (size_t)k)];
        };

        double c00 = at(i0,   j0,   k0  ) * (1.0 - di) + at(i0+1, j0,   k0  ) * di;
        double c10 = at(i0,   j0+1, k0  ) * (1.0 - di) + at(i0+1, j0+1, k0  ) * di;
        double c01 = at(i0,   j0,   k0+1) * (1.0 - di) + at(i0+1, j0,   k0+1) * di;
        double c11 = at(i0,   j0+1, k0+1) * (1.0 - di) + at(i0+1, j0+1, k0+1) * di;

        double c0 = c00 * (1.0 - dj) + c10 * dj;
        double c1 = c01 * (1.0 - dj) + c11 * dj;

        out = (float)(c0 * (1.0 - dk) + c1 * dk);
        return true;
    }
};

// ---------------------------------------------------------------------------
// FORCE TERM 1: intensity-target localization.
//
// For each vertex, samples the volume along the vertex's *current* normal at
// offsets -max_thickness, -max_thickness + step, ..., +max_thickness, and
// picks the offset whose sampled intensity is closest to `target`. That
// offset is the displacement (along the normal) that would place the vertex
// where the local intensity profile best matches the target; it is scaled by
// `l_intensity` and added directly to the gradient.
//
// This is a deliberate simplification of the search procedure described in
// the literature, which additionally reasons about profile monotonicity and
// derivatives, and adapts the target per vertex from segmentation-derived
// gray/white/'CSF' statistics. "Closest match to a single fixed target"
// keeps the same underlying idea, pull the surface to the depth along its
// normal where the tissue-intensity boundary is, with one tunable number per
// surface rather than a per-vertex statistical model.
// ---------------------------------------------------------------------------
void intensityTargetTerm(MrisMesh &m, const IntensityVolume &vol, float target,
                         float max_thickness, float step, float l_intensity)
{
    if (l_intensity == 0.0f || step <= 0.0f) return;
    int nsteps = (int)std::floor(max_thickness / step);

    for (int vi = 0; vi < m.nv; vi++) {
        float px = m.x[vi], py = m.y[vi], pz = m.z[vi];
        float vnx = m.nx[vi], vny = m.ny[vi], vnz = m.nz[vi];

        float best_t = 0.0f;
        float best_diff = -1.0f;
        bool found = false;

        for (int s = -nsteps; s <= nsteps; s++) {
            float t = s * step;
            float val;
            if (!vol.sample(px + t*vnx, py + t*vny, pz + t*vnz, val)) continue;
            float diff = std::fabs(val - target);
            if (best_diff < 0.0f || diff < best_diff) { best_diff = diff; best_t = t; found = true; }
        }
        if (!found) continue;

        m.dx[vi] += l_intensity * best_t * vnx;
        m.dy[vi] += l_intensity * best_t * vny;
        m.dz[vi] += l_intensity * best_t * vnz;
    }
}

// ---------------------------------------------------------------------------
// FORCE TERM 2: 1-ring Laplacian smoothness.
//
// Pulls each vertex toward the centroid of its 1-ring neighbors: a plain
// mean-curvature-flow-style spring that keeps the mesh regular while the
// intensity term (which acts independently per vertex from noisy image data)
// pulls individual vertices toward the tissue boundary.
// ---------------------------------------------------------------------------
void surfaceSmoothnessTerm(MrisMesh &m, float l_spring)
{
    if (l_spring == 0.0f) return;

    for (int vi = 0; vi < m.nv; vi++) {
        const auto &nb = m.nbrs1[vi];
        int nn = (int)nb.size();
        if (nn == 0) continue;

        float sx = 0.0f, sy = 0.0f, sz = 0.0f;
        for (int nj : nb) { sx += m.x[nj]; sy += m.y[nj]; sz += m.z[nj]; }
        float inv = 1.0f / nn;

        m.dx[vi] += l_spring * (sx * inv - m.x[vi]);
        m.dy[vi] += l_spring * (sy * inv - m.y[vi]);
        m.dz[vi] += l_spring * (sz * inv - m.z[vi]);
    }
}

// ---------------------------------------------------------------------------
// Per-surface deformation pass. `niterations` rounds of:
//   clear gradient -> intensity-target term -> average gradient (n_averages
//   passes of 1-ring smoothing, the same gradient-averaging building block
//   mris_inflate / mris_sphere use) -> smoothness term (added after
//   averaging, acts locally, mirroring how mris_inflate adds its spring term)
//   -> momentum step -> refresh normals (the only per-vertex quantity this
//   loop's force terms depend on; areas/distances are not used).
// ---------------------------------------------------------------------------
void deformTowardIntensityTarget(MrisMesh &m, const IntensityVolume &vol, float target,
                                 float max_thickness, float step_size,
                                 int n_averages, int niterations,
                                 float l_intensity, float l_spring,
                                 float momentum, float dt,
                                 const char *label, bool verbose)
{
    for (int n = 0; n < niterations; n++) {
        Rcpp::checkUserInterrupt();

        std::fill(m.dx.begin(), m.dx.end(), 0.0f);
        std::fill(m.dy.begin(), m.dy.end(), 0.0f);
        std::fill(m.dz.begin(), m.dz.end(), 0.0f);

        intensityTargetTerm(m, vol, target, max_thickness, step_size, l_intensity);
        mrisAverageField(m, m.dx, m.dy, m.dz, n_averages);
        surfaceSmoothnessTerm(m, l_spring);

        mrisMomentumStep(m, momentum, dt, MAX_SURFACE_STEP_MM);
        mrisComputeNormals(m);

        if (verbose && ((n + 1) % 5 == 0)) {
            Rprintf("  %s surface: iter %3d/%d\n", label, n + 1, niterations);
        } else {
            Rcpp::checkUserInterrupt();
        }
    }
}

// ===========================================================================
// mris_curvature helpers
// ===========================================================================

// ---------------------------------------------------------------------------
// Solve a 3x3 linear system A x = b via Cramer's rule. Returns false (leaving
// `x` untouched) when A is (numerically) singular, e.g. a vertex whose
// neighborhood is degenerate (collinear / too few neighbors).
// ---------------------------------------------------------------------------
bool solveLinear3x3(const double A[3][3], const double b[3], double x[3])
{
    auto det3 = [](double m00, double m01, double m02,
                   double m10, double m11, double m12,
                   double m20, double m21, double m22) -> double {
        return m00*(m11*m22 - m12*m21)
             - m01*(m10*m22 - m12*m20)
             + m02*(m10*m21 - m11*m20);
    };

    double det = det3(A[0][0], A[0][1], A[0][2],
                      A[1][0], A[1][1], A[1][2],
                      A[2][0], A[2][1], A[2][2]);
    if (std::fabs(det) < 1e-12) return false;

    // Cramer's rule: x[col] = det(A with column `col` replaced by b) / det(A)
    x[0] = det3(b[0],    A[0][1], A[0][2],
                b[1],    A[1][1], A[1][2],
                b[2],    A[2][1], A[2][2]) / det;
    x[1] = det3(A[0][0], b[0],    A[0][2],
                A[1][0], b[1],    A[1][2],
                A[2][0], b[2],    A[2][2]) / det;
    x[2] = det3(A[0][0], A[0][1], b[0],
                A[1][0], A[1][1], b[1],
                A[2][0], A[2][1], b[2]) / det;
    return true;
}

// ---------------------------------------------------------------------------
// Per-vertex curvature via local quadratic (osculating-paraboloid) fitting,
// the discrete-differential-geometry analogue of the literature's
// second-fundamental-form estimation:
//
//   1. Build an orthonormal tangent frame (e1, e2, n) at the vertex, with `n`
//      its (already-computed) unit normal.
//   2. Express each 2-ring neighbor's offset (p - v) in this frame as
//      (u, w, h): `u`, `w` are tangential coordinates, `h` the height above
//      the tangent plane along `n`.
//   3. Least-squares fit the osculating paraboloid h = a*u^2 + b*u*w + c*w^2
//      (no linear/constant terms: by construction the surface passes through
//      the origin with zero gradient there, since `n` is the fitted normal).
//   4. Because (e1, e2) is orthonormal, the second fundamental form is
//      [[2a, b], [b, 2c]] and the first fundamental form is the identity, so
//      its eigenvalues ARE the principal curvatures:
//        mean curvature     H = (k1 + k2)/2 = a + c
//        Gaussian curvature K =  k1 * k2    = 4*a*c - b^2
//        k1, k2 = H +/- sqrt(max(H^2 - K, 0))
//
// Sign convention follows the orientation of the per-vertex outward normal
// `n`: neighbors that lie on the side `n` points away from (h < 0, the usual
// case for an outward normal on a locally convex / "gyrus-like" patch) yield
// negative curvature; neighbors on the side `n` points toward (h > 0, a
// locally concave / "sulcus-like" patch) yield positive curvature. A sphere
// of radius r with outward normals is uniformly curved with H = -1/r,
// K = 1/r^2.
//
// Vertices whose neighborhood fit is degenerate (fewer than 3 neighbors, or
// a singular normal-equations matrix, e.g. collinear neighbors) are left at
// zero.
// ---------------------------------------------------------------------------
void computeVertexCurvatures(const MrisMesh &m,
                             std::vector<float> &mean_curv, std::vector<float> &gauss_curv,
                             std::vector<float> &k1, std::vector<float> &k2)
{
    mean_curv.assign(m.nv, 0.0f);
    gauss_curv.assign(m.nv, 0.0f);
    k1.assign(m.nv, 0.0f);
    k2.assign(m.nv, 0.0f);

    for (int vi = 0; vi < m.nv; vi++) {
        float nxv = m.nx[vi], nyv = m.ny[vi], nzv = m.nz[vi];

        // Pick the coordinate axis least aligned with the normal as the seed
        // for building an orthonormal tangent basis (avoids a near-zero cross
        // product).
        float ax, ay, az;
        float anx = std::fabs(nxv), any = std::fabs(nyv), anz = std::fabs(nzv);
        if (anx <= any && anx <= anz)      { ax = 1.0f; ay = 0.0f; az = 0.0f; }
        else if (any <= anz)               { ax = 0.0f; ay = 1.0f; az = 0.0f; }
        else                               { ax = 0.0f; ay = 0.0f; az = 1.0f; }

        float e1x = nyv*az - nzv*ay, e1y = nzv*ax - nxv*az, e1z = nxv*ay - nyv*ax;
        float e1l = std::sqrt(e1x*e1x + e1y*e1y + e1z*e1z);
        e1x /= e1l; e1y /= e1l; e1z /= e1l;

        float e2x = nyv*e1z - nzv*e1y, e2y = nzv*e1x - nxv*e1z, e2z = nxv*e1y - nyv*e1x;

        const auto &nb = m.nbrs2[vi];
        int nn = (int)nb.size();
        if (nn < 3) continue;

        double AtA[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        double Atb[3] = {0.0, 0.0, 0.0};

        for (int ni = 0; ni < nn; ni++) {
            int nj = nb[ni];
            float ddx = m.x[nj]-m.x[vi], ddy = m.y[nj]-m.y[vi], ddz = m.z[nj]-m.z[vi];

            double u = (double)(ddx*e1x + ddy*e1y + ddz*e1z);
            double w = (double)(ddx*e2x + ddy*e2y + ddz*e2z);
            double h = (double)(ddx*nxv + ddy*nyv + ddz*nzv);

            double f[3] = { u*u, u*w, w*w };
            for (int r = 0; r < 3; r++) {
                Atb[r] += f[r] * h;
                for (int c = 0; c < 3; c++) AtA[r][c] += f[r] * f[c];
            }
        }

        double abc[3];
        if (!solveLinear3x3(AtA, Atb, abc)) continue;

        double a = abc[0], b = abc[1], c = abc[2];
        double Hv = a + c;
        double Kv = 4.0*a*c - b*b;
        double disc = Hv*Hv - Kv;
        if (disc < 0.0) disc = 0.0;
        double sq = std::sqrt(disc);

        mean_curv[vi]  = (float)Hv;
        gauss_curv[vi] = (float)Kv;
        k1[vi]         = (float)(Hv + sq);
        k2[vi]         = (float)(Hv - sq);
    }
}

// ===========================================================================
// mris_remesh helpers (Botsch & Kobbelt 2004 isotropic remeshing)
// ===========================================================================

// ---------------------------------------------------------------------------
// Canonical 64-bit key for an undirected edge (a, b), with a < b.
// Same encoding as mrisCheckClosedManifold's edge_count map.
// ---------------------------------------------------------------------------
inline int64_t remeshEdgeKey(int a, int b)
{
    if (a > b) { int t = a; a = b; b = t; }
    return (int64_t)a << 32 | (int64_t)(uint32_t)b;
}

// ---------------------------------------------------------------------------
// Edge split: for every edge longer than `max_len`, insert a midpoint vertex
// and replace the face(s) sharing that edge with two sub-triangles each.
//
// 1-edge split (one long edge per face):  1 face -> 2
// 2-edge split (two long edges):           1 face -> 3
// 3-edge split (all three long edges):     1 face -> 4  (standard 1-to-4)
//
// Per-vertex ancillary arrays (nx,ny,nz,dx,dy,dz,odx,ody,odz) are zero-
// initialised for new midpoint vertices; mrisComputeNormals resets them
// immediately after this call.
// ---------------------------------------------------------------------------
void mrisRemeshSplitEdges(MrisMesh &m, float max_len)
{
    const float max_len2 = max_len * max_len;

    // Map edge -> new midpoint vertex index
    std::unordered_map<int64_t, int> midVert;
    midVert.reserve((size_t)m.nf * 3 / 2);

    // Pass 1: find long edges, allocate midpoint vertices
    auto edgeMid = [&](int a, int b) -> int {
        int64_t key = remeshEdgeKey(a, b);
        auto it = midVert.find(key);
        if (it != midVert.end()) return it->second;

        float dx = m.x[b]-m.x[a], dy = m.y[b]-m.y[a], dz = m.z[b]-m.z[a];
        if (dx*dx + dy*dy + dz*dz <= max_len2) return -1;

        int nv = (int)m.x.size();
        m.x.push_back(0.5f*(m.x[a]+m.x[b]));
        m.y.push_back(0.5f*(m.y[a]+m.y[b]));
        m.z.push_back(0.5f*(m.z[a]+m.z[b]));
        m.nx.push_back(0.0f); m.ny.push_back(0.0f); m.nz.push_back(0.0f);
        m.dx.push_back(0.0f); m.dy.push_back(0.0f); m.dz.push_back(0.0f);
        m.odx.push_back(0.0f); m.ody.push_back(0.0f); m.odz.push_back(0.0f);
        midVert[key] = nv;
        return nv;
    };

    int orig_nf = m.nf;
    for (int fi = 0; fi < orig_nf; fi++) {
        edgeMid(m.f0[fi], m.f1[fi]);
        edgeMid(m.f1[fi], m.f2[fi]);
        edgeMid(m.f2[fi], m.f0[fi]);
    }

    // Pass 2: split faces
    // Reserve extra space (worst case 4x faces)
    m.f0.reserve(m.f0.size() * 4);
    m.f1.reserve(m.f1.size() * 4);
    m.f2.reserve(m.f2.size() * 4);

    for (int fi = 0; fi < orig_nf; fi++) {
        int A = m.f0[fi], B = m.f1[fi], C = m.f2[fi];
        int mAB = midVert.count(remeshEdgeKey(A,B)) ? midVert[remeshEdgeKey(A,B)] : -1;
        int mBC = midVert.count(remeshEdgeKey(B,C)) ? midVert[remeshEdgeKey(B,C)] : -1;
        int mCA = midVert.count(remeshEdgeKey(C,A)) ? midVert[remeshEdgeKey(C,A)] : -1;

        bool sAB = (mAB >= 0), sBC = (mBC >= 0), sCA = (mCA >= 0);
        int n_split = (int)sAB + (int)sBC + (int)sCA;

        if (n_split == 0) continue;

        if (n_split == 3) {
            // Standard 1-to-4 split
            m.f0[fi] = A;   m.f1[fi] = mAB; m.f2[fi] = mCA;
            m.f0.push_back(mAB); m.f1.push_back(B);   m.f2.push_back(mBC);
            m.f0.push_back(mCA); m.f1.push_back(mBC); m.f2.push_back(C);
            m.f0.push_back(mAB); m.f1.push_back(mBC); m.f2.push_back(mCA);
        } else if (n_split == 1) {
            if (sAB) {
                m.f0[fi] = A;   m.f1[fi] = mAB; m.f2[fi] = C;
                m.f0.push_back(mAB); m.f1.push_back(B); m.f2.push_back(C);
            } else if (sBC) {
                m.f0[fi] = B;   m.f1[fi] = mBC; m.f2[fi] = A;
                m.f0.push_back(mBC); m.f1.push_back(C); m.f2.push_back(A);
            } else {
                m.f0[fi] = C;   m.f1[fi] = mCA; m.f2[fi] = B;
                m.f0.push_back(mCA); m.f1.push_back(A); m.f2.push_back(B);
            }
        } else {
            // n_split == 2: split the two cut edges; the third is the longest sub-edge
            // Arrange so the single un-split edge is CA when we enter the logic:
            //   rotate ABC so the un-split edge is always between C and A.
            int vA = A, vB = B, vC = C;
            int vmAB = mAB, vmBC = mBC;
            if (!sAB) {
                // un-split: AB -> rotate: A=B, B=C, C=A
                vA = B; vB = C; vC = A; vmAB = mBC; vmBC = mCA;
            } else if (!sBC) {
                // un-split: BC -> rotate: A=C, B=A, C=B
                vA = C; vB = A; vC = B; vmAB = mCA; vmBC = mAB;
            }
            // Now vA-vB and vB-vC are split, vC-vA is not.
            // Split into 3 triangles:
            //   (vA, vmAB, vC), (vmAB, vmBC, vC), (vmAB, vB, vmBC)
            m.f0[fi] = vA;   m.f1[fi] = vmAB; m.f2[fi] = vC;
            m.f0.push_back(vmAB); m.f1.push_back(vmBC); m.f2.push_back(vC);
            m.f0.push_back(vmAB); m.f1.push_back(vB);   m.f2.push_back(vmBC);
        }
    }

    m.nv = (int)m.x.size();
    m.nf = (int)m.f0.size();
}

// ---------------------------------------------------------------------------
// Edge collapse: for every edge shorter than `min_len`, if collapsing it
// (moving vertex a to vertex b's position, then removing a) would not flip
// any surrounding face normal, do the collapse. After all collapses in a
// single pass, compact the vertex and face arrays.
// ---------------------------------------------------------------------------
void mrisRemeshCollapseEdges(MrisMesh &m, float min_len)
{
    const float min_len2 = min_len * min_len;

    // Collect vertex-to-face adjacency (1-ring faces per vertex)
    std::vector<std::vector<int>> v2f(m.nv);
    for (int fi = 0; fi < m.nf; fi++) {
        v2f[m.f0[fi]].push_back(fi);
        v2f[m.f1[fi]].push_back(fi);
        v2f[m.f2[fi]].push_back(fi);
    }

    std::vector<bool> collapsed(m.nv, false);
    std::vector<bool> face_del(m.nf, false);

    // Helper: cross product z-component sign of triangle (p0,p1,p2) projected
    // onto the plane with normal n; returns true if the winding is consistent
    // with n (dot(cross, n) > 0).
    auto normalOk = [&](int v0, int v1, int v2, int ref_normal_v) -> bool {
        float ax = m.x[v1]-m.x[v0], ay = m.y[v1]-m.y[v0], az = m.z[v1]-m.z[v0];
        float bx = m.x[v2]-m.x[v0], by = m.y[v2]-m.y[v0], bz = m.z[v2]-m.z[v0];
        float cx = ay*bz - az*by, cy = az*bx - ax*bz, cz = ax*by - ay*bx;
        return cx*m.nx[ref_normal_v] + cy*m.ny[ref_normal_v] + cz*m.nz[ref_normal_v] > 0.0f;
    };

    // Enumerate unique edges (short ones)
    std::unordered_map<int64_t, std::pair<int,int>> edges;
    edges.reserve((size_t)m.nf * 3 / 2);
    for (int fi = 0; fi < m.nf; fi++) {
        int vs[3] = { m.f0[fi], m.f1[fi], m.f2[fi] };
        for (int e = 0; e < 3; e++) {
            int a = vs[e], b = vs[(e+1)%3];
            int64_t key = remeshEdgeKey(a, b);
            edges.emplace(key, std::make_pair(a, b));
        }
    }

    for (auto &kv : edges) {
        int a = kv.second.first, b = kv.second.second;
        if (collapsed[a] || collapsed[b]) continue;

        float ddx = m.x[b]-m.x[a], ddy = m.y[b]-m.y[a], ddz = m.z[b]-m.z[a];
        if (ddx*ddx + ddy*ddy + ddz*ddz > min_len2) continue;

        // Target position: midpoint
        float mx = 0.5f*(m.x[a]+m.x[b]);
        float my = 0.5f*(m.y[a]+m.y[b]);
        float mz = 0.5f*(m.z[a]+m.z[b]);

        // Tentatively move b to midpoint; check no face flips around a or b
        float old_bx = m.x[b], old_by = m.y[b], old_bz = m.z[b];
        m.x[b] = mx; m.y[b] = my; m.z[b] = mz;

        bool valid = true;
        for (int fi : v2f[b]) {
            if (face_del[fi]) continue;
            int v0 = m.f0[fi], v1 = m.f1[fi], v2 = m.f2[fi];
            // Skip faces that share edge (a,b) — they will be deleted
            if ((v0==a||v1==a||v2==a) && (v0==b||v1==b||v2==b)) continue;
            // Skip degenerate
            if (v0==v1 || v1==v2 || v0==v2) continue;
            if (!normalOk(v0, v1, v2, b)) { valid = false; break; }
        }
        if (valid) {
            for (int fi : v2f[a]) {
                if (face_del[fi]) continue;
                int v0 = m.f0[fi], v1 = m.f1[fi], v2 = m.f2[fi];
                if ((v0==a||v1==a||v2==a) && (v0==b||v1==b||v2==b)) continue;
                if (v0==v1 || v1==v2 || v0==v2) continue;
                // Temporarily remap a->b in this face for the normal check
                int r0 = (v0==a)?b:v0, r1 = (v1==a)?b:v1, r2 = (v2==a)?b:v2;
                if (r0==r1 || r1==r2 || r0==r2) continue;
                if (!normalOk(r0, r1, r2, b)) { valid = false; break; }
            }
        }

        if (!valid) {
            m.x[b] = old_bx; m.y[b] = old_by; m.z[b] = old_bz;
            continue;
        }

        // Commit collapse: a is removed, b takes the midpoint position
        collapsed[a] = true;
        for (int fi : v2f[a]) {
            if (face_del[fi]) continue;
            int &r0 = m.f0[fi], &r1 = m.f1[fi], &r2 = m.f2[fi];
            if (r0==a) r0 = b;
            if (r1==a) r1 = b;
            if (r2==a) r2 = b;
            if (r0==r1 || r1==r2 || r0==r2) face_del[fi] = true;
            else v2f[b].push_back(fi);
        }
        // Mark faces that were already on edge (a,b) as deleted
        for (int fi : v2f[b]) {
            if (!face_del[fi]) {
                int v0=m.f0[fi], v1=m.f1[fi], v2=m.f2[fi];
                if (v0==v1 || v1==v2 || v0==v2) face_del[fi] = true;
            }
        }
    }

    // Compact: build alive vertex set and remap
    std::vector<int> remap(m.nv, -1);
    int new_nv = 0;
    for (int vi = 0; vi < m.nv; vi++) {
        if (!collapsed[vi]) remap[vi] = new_nv++;
    }

    std::vector<float> nx(new_nv), ny(new_nv), nz(new_nv);
    std::vector<float> xx(new_nv), xy(new_nv), xz(new_nv);
    std::vector<float> ddx(new_nv,0), ddy(new_nv,0), ddz(new_nv,0);
    std::vector<float> odx(new_nv,0), ody(new_nv,0), odz(new_nv,0);
    for (int vi = 0; vi < m.nv; vi++) {
        int ri = remap[vi]; if (ri < 0) continue;
        xx[ri]=m.x[vi]; xy[ri]=m.y[vi]; xz[ri]=m.z[vi];
        nx[ri]=m.nx[vi]; ny[ri]=m.ny[vi]; nz[ri]=m.nz[vi];
    }
    m.x=xx; m.y=xy; m.z=xz;
    m.nx=nx; m.ny=ny; m.nz=nz;
    m.dx=ddx; m.dy=ddy; m.dz=ddz;
    m.odx=odx; m.ody=ody; m.odz=odz;
    m.nv = new_nv;

    std::vector<int> nf0, nf1, nf2;
    nf0.reserve(m.nf); nf1.reserve(m.nf); nf2.reserve(m.nf);
    for (int fi = 0; fi < m.nf; fi++) {
        if (face_del[fi]) continue;
        int r0 = remap[m.f0[fi]], r1 = remap[m.f1[fi]], r2 = remap[m.f2[fi]];
        if (r0<0||r1<0||r2<0) continue;
        nf0.push_back(r0); nf1.push_back(r1); nf2.push_back(r2);
    }
    m.f0=nf0; m.f1=nf1; m.f2=nf2;
    m.nf = (int)m.f0.size();
}

// ---------------------------------------------------------------------------
// Tangential smooth: `n_iters` passes of uniform-weight 1-ring Laplacian,
// with each displacement projected onto the vertex tangent plane before
// application (preserves surface shape). Requires m.nbrs1 populated.
// ---------------------------------------------------------------------------
void mrisRemeshTangentialSmooth(MrisMesh &m, int n_iters, float damping)
{
    for (int iter = 0; iter < n_iters; iter++) {
        for (int vi = 0; vi < m.nv; vi++) {
            const auto &nb = m.nbrs1[vi];
            int nn = (int)nb.size();
            if (nn == 0) continue;

            float cx = 0, cy = 0, cz = 0;
            for (int nj : nb) { cx += m.x[nj]; cy += m.y[nj]; cz += m.z[nj]; }
            cx /= nn; cy /= nn; cz /= nn;

            float dispx = cx - m.x[vi], dispy = cy - m.y[vi], dispz = cz - m.z[vi];
            // Remove normal component (project onto tangent plane)
            float nc = dispx*m.nx[vi] + dispy*m.ny[vi] + dispz*m.nz[vi];
            dispx -= nc*m.nx[vi]; dispy -= nc*m.ny[vi]; dispz -= nc*m.nz[vi];

            m.x[vi] += damping * dispx;
            m.y[vi] += damping * dispy;
            m.z[vi] += damping * dispz;
        }
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

    // mrisCheckClosedManifold(m, "mris_smooth");

    if (verbose) {
        Rprintf("mris_smooth: building adjacency for %d vertices, %d faces\n", m.nv, m.nf);
    } else {
        Rcpp::checkUserInterrupt();
    }
    mrisBuildAdjacency(m, /* two_ring = */ false);

    mrisComputeNormals(m);
    mrisComputeArea(m);
    m.orig_area = m.total_area;   // reference area, used to rescale back if requested

    for (int pass = 0; pass < npasses; pass++) {
        Rcpp::checkUserInterrupt();

        if (verbose) {
            Rprintf("mris_smooth: pass %d/%d: averaging vertex positions over %d iterations\n",
                    pass + 1, npasses, niterations);
        }

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

    if (verbose) {
        Rprintf("mris_smooth: done\n");
    }

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
    // mrisCheckClosedManifold(m, "mris_inflate");

    std::vector<float> sulc(m.nv, 0.0f);   // per-vertex sulcal depth accumulator

    // --- Pre-processing ---

    if (verbose) {
        Rprintf("mris_inflate: building adjacency for %d vertices, %d faces\n", m.nv, m.nf);
    }
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
    if (verbose) {
        Rprintf("mris_inflate: computing original distances\n");
    }
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

    if (verbose) {
        Rprintf("mris_inflate: done (%d total iterations)\n", total_iters);
    }

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

    // mrisCheckClosedManifold(m, "mris_sphere");

    if (verbose) {
        Rprintf("mris_sphere: building adjacency for %d vertices, %d faces\n", m.nv, m.nf);
    }
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
    if (verbose) {
        Rprintf("mris_sphere: projecting onto sphere of radius %.2f\n", radius);
    }
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

    if (verbose) {
        Rprintf("mris_sphere: done (%d iterations)\n", niterations);
    }

    return mrisPackMeshOutput(m);
}

// ===========================================================================
// mris_make_surfaces
//
// Intensity-driven surface placement: starting from a single closed surface
// (typically a smoothed estimate of the white-matter boundary) and a
// co-registered intensity volume, deforms it in two passes to localize the
// white/gray and gray/'CSF' tissue-intensity boundaries, producing a "white"
// and a "pial" surface.
//
// This is a *reduced* implementation, not a faithful reproduction of the
// literature's procedure: that drives a multi-resolution optimization over
// roughly seven weighted energy terms (intensity, intensity gradient,
// smoothness, self-intersection repulsion, curvature, ...), using per-vertex
// gray/white/'CSF' intensity statistics derived from a prior segmentation.
// Reproducing that faithfully is out of scope for this package. Instead this
// keeps the two dominant ideas:
//   1. an intensity-target localization force that searches the volume along
//      each vertex's normal for the offset whose intensity best matches a
//      single caller-supplied target value (see intensityTargetTerm), and
//   2. a 1-ring Laplacian smoothness force that keeps the mesh regular
//      (see surfaceSmoothnessTerm),
// integrated with the same gradient-averaging / momentum-integration
// machinery that mris_inflate and mris_sphere use. A single target intensity
// is supplied per surface (`white_intensity`, `pial_intensity`) in place of
// the segmentation-derived per-vertex statistics, e.g. the midpoint between
// the typical white-matter and gray-matter intensities for the white surface,
// and between gray-matter and 'CSF' for the pial surface.
//
// Reference: Dale AM, Fischl B, Sereno MI (1999). Cortical surface-based
// analysis I: Segmentation and surface reconstruction. NeuroImage, 9(2),
// 179-194.
// ===========================================================================

// ---------------------------------------------------------------------------
// Main exported function.
//
// Pipeline:
//   1. Build 1-ring adjacency and normals for the input surface (no scaling
//      or centering: the surface and volume must already share the same
//      physical 'RAS' coordinate space).
//   2. White pass: deform the input surface toward `white_intensity`.
//   3. Pial pass: continue from the resulting white surface, with the
//      momentum velocity reset to rest, deforming toward `pial_intensity`
//      along the (recomputed-each-iteration) current normal; since the pial
//      boundary normally lies further outward from the white surface along
//      the local normal, this is the natural continuation of the same search.
//
// `volume` is a 3-D numeric array; `ras2ijk` the 4x4 affine mapping a
// surface-space ('RAS') point to fractional volume 'IJK' indices (the
// inverse of the volume's 'IJK'-to-'RAS' transform). Both are prepared by the
// R wrapper.
//
// Returns a list with `white` and `pial` entries, each itself a list with
// `vb` (3 x nv positions) and `normals` (3 x nv normals), matching the
// per-surface output shape of mris_inflate / mris_sphere.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mrisMakeSurfaces(
    SEXP   vb_,
    SEXP   it_,
    SEXP   volume_,
    SEXP   ras2ijk_,
    double white_intensity,
    double pial_intensity,
    double max_thickness = 5.0,
    double step_size     = 0.4,
    int    n_averages    = 4,
    int    niterations   = 10,
    double l_intensity   = 1.0,
    double l_spring      = 0.5,
    double momentum      = 0.9,
    double dt            = 0.5,
    bool   verbose       = false
) {
    MrisMesh m;
    mrisLoadMesh(m, vb_, it_);
    // mrisCheckClosedManifold(m, "mris_make_surfaces");

    if (verbose) {
        Rprintf("mris_make_surfaces: building adjacency for %d vertices, %d faces\n",
                m.nv, m.nf);
    }
    mrisBuildAdjacency(m, /* two_ring = */ false);
    mrisComputeNormals(m);

    // --- Intensity volume sampler ---
    Rcpp::NumericVector vol(volume_);
    Rcpp::IntegerVector vol_dim = vol.attr("dim");
    if (vol_dim.size() != 3) {
        Rcpp::stop("mris_make_surfaces: `volume` must be a 3-dimensional array");
    }
    Rcpp::NumericMatrix ras2ijk(ras2ijk_);
    if (ras2ijk.nrow() != 4 || ras2ijk.ncol() != 4) {
        Rcpp::stop("mris_make_surfaces: `ras2ijk` must be a 4x4 matrix");
    }

    IntensityVolume vol_sampler;
    vol_sampler.data = vol.begin();
    vol_sampler.nx   = vol_dim[0];
    vol_sampler.ny   = vol_dim[1];
    vol_sampler.nz   = vol_dim[2];
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 4; c++) {
            vol_sampler.r2i[r*4 + c] = (float)ras2ijk(r, c);
        }
    }

    float max_thickness_f = (float)max_thickness;
    float step_size_f     = (float)step_size;
    float l_intensity_f   = (float)l_intensity;
    float l_spring_f      = (float)l_spring;
    float momentum_f      = (float)momentum;
    float dt_f            = (float)dt;

    // --- Pass 1: white surface ---
    if (verbose) {
      Rprintf("mris_make_surfaces: localizing white surface (target intensity = %.2f)\n",
              white_intensity);
    }
    deformTowardIntensityTarget(m, vol_sampler, (float)white_intensity,
                                max_thickness_f, step_size_f,
                                n_averages, niterations,
                                l_intensity_f, l_spring_f, momentum_f, dt_f,
                                "white", verbose);
    Rcpp::List white_out = mrisPackMeshOutput(m);

    // --- Pial pass: continue from the white surface; reset momentum to rest ---
    if (verbose) {
      Rprintf("mris_make_surfaces: localizing pial surface (target intensity = %.2f)\n",
              pial_intensity);
    }
    std::fill(m.odx.begin(), m.odx.end(), 0.0f);
    std::fill(m.ody.begin(), m.ody.end(), 0.0f);
    std::fill(m.odz.begin(), m.odz.end(), 0.0f);

    deformTowardIntensityTarget(m, vol_sampler, (float)pial_intensity,
                                max_thickness_f, step_size_f,
                                n_averages, niterations,
                                l_intensity_f, l_spring_f, momentum_f, dt_f,
                                "pial", verbose);
    Rcpp::List pial_out = mrisPackMeshOutput(m);

    if (verbose) {
        Rprintf("mris_make_surfaces: done\n");
    }

    return Rcpp::List::create(
        Rcpp::Named("white") = white_out,
        Rcpp::Named("pial")  = pial_out
    );
}

// ===========================================================================
// mris_curvature
// ===========================================================================
//
// Per-vertex curvature estimates for a closed surface mesh, via local
// quadratic (osculating-paraboloid) fitting over each vertex's 2-ring
// neighborhood; see computeVertexCurvatures above for the full algorithm and
// sign-convention notes.
//
// Returns four numeric vectors of length `nv`: `mean` (mean curvature `H`),
// `gaussian` (Gaussian curvature `K`), and the two principal curvatures
// `k1`, `k2` (with `k1 >= k2` and `H = (k1+k2)/2`, `K = k1*k2`).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mrisCurvature(
    SEXP vb_,
    SEXP it_,
    bool verbose = false
) {
    MrisMesh m;
    mrisLoadMesh(m, vb_, it_);
    // mrisCheckClosedManifold(m, "mris_curvature");

    if (verbose) {
      Rprintf("mris_curvature: building adjacency for %d vertices, %d faces\n",
              m.nv, m.nf);
    }
    mrisBuildAdjacency(m, /* two_ring = */ true);
    mrisComputeNormals(m);

    if (verbose) {
      Rprintf("mris_curvature: fitting local quadratics over 2-ring neighborhoods\n");
    }
    std::vector<float> mean_curv, gauss_curv, k1, k2;
    computeVertexCurvatures(m, mean_curv, gauss_curv, k1, k2);

    if (verbose) {
      Rprintf("mris_curvature: done\n");
    }

    return Rcpp::List::create(
        Rcpp::Named("mean")     = Rcpp::wrap(mean_curv),
        Rcpp::Named("gaussian") = Rcpp::wrap(gauss_curv),
        Rcpp::Named("k1")       = Rcpp::wrap(k1),
        Rcpp::Named("k2")       = Rcpp::wrap(k2)
    );
}

// ===========================================================================
// mris_remesh
// ===========================================================================
//
// Isotropic remeshing following the Botsch & Kobbelt 2004 algorithm: each
// iteration applies edge splitting (for edges > 4/3 * target), edge collapse
// (for edges < 4/5 * target), and tangential smoothing (uniform-weight
// Laplacian projected onto the tangent plane, preserving surface shape).
//
// Unlike vcg_uniform_remesh (volumetric resampling), this operates purely on
// the mesh surface and preserves the input topology class.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mrisRemesh(
    SEXP  vb_,
    SEXP  it_,
    float target_edge_len,
    int   niterations = 5,
    int   n_smooth    = 2,
    float damping     = 0.99,
    bool  verbose     = false
) {
    MrisMesh m;
    mrisLoadMesh(m, vb_, it_);
    mrisCheckClosedManifold(m, "mris_remesh");

    mrisComputeNormals(m);

    const float max_len = (4.0f / 3.0f) * target_edge_len;
    const float min_len = (4.0f / 5.0f) * target_edge_len;

    if (verbose) {
        Rprintf("mris_remesh: target=%.4f  split>%.4f  collapse<%.4f\n",
                target_edge_len, max_len, min_len);
    }

    for (int iter = 0; iter < niterations; iter++) {
        mrisRemeshSplitEdges(m, max_len);
        mrisBuildAdjacency(m, /* two_ring = */ false);
        mrisComputeNormals(m);

        mrisRemeshCollapseEdges(m, min_len);
        mrisBuildAdjacency(m, /* two_ring = */ false);
        mrisComputeNormals(m);

        mrisRemeshTangentialSmooth(m, n_smooth, damping);
        mrisComputeNormals(m);

        if (verbose) {
            Rprintf("  iter %d: nv=%d  nf=%d\n", iter + 1, m.nv, m.nf);
        }
        Rcpp::checkUserInterrupt();
    }

    // Pack positions and normals
    Rcpp::List mesh_out = mrisPackMeshOutput(m);

    // Pack face indices (1-based for R's mesh3d convention)
    Rcpp::IntegerMatrix out_it(3, m.nf);
    for (int fi = 0; fi < m.nf; fi++) {
        out_it(0, fi) = m.f0[fi] + 1;
        out_it(1, fi) = m.f1[fi] + 1;
        out_it(2, fi) = m.f2[fi] + 1;
    }

    return Rcpp::List::create(
        Rcpp::Named("vb")      = mesh_out["vb"],
        Rcpp::Named("normals") = mesh_out["normals"],
        Rcpp::Named("it")      = out_it
    );
}
