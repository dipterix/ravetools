// mrisCommon.cpp
//
// Shared implementation for the cortical-surface algorithm ports.
// See mrisCommon.h for the rationale and references.

#include "mrisCommon.h"
#include <algorithm>
#include <unordered_map>
#include <cstdint>

namespace ravetools {

// ---------------------------------------------------------------------------
void mrisLoadMesh(MrisMesh &m, SEXP vb_, SEXP it_)
{
    Rcpp::NumericMatrix vb(vb_);
    Rcpp::IntegerMatrix it(it_);

    m.nv = vb.ncol();
    m.nf = it.ncol();

    if (m.nv < 4 || m.nf < 1) {
        Rcpp::stop("mris: mesh has too few vertices or faces");
    }

    m.x.resize(m.nv); m.y.resize(m.nv); m.z.resize(m.nv);
    for (int vi = 0; vi < m.nv; vi++) {
        m.x[vi] = (float)vb(0, vi);
        m.y[vi] = (float)vb(1, vi);
        m.z[vi] = (float)vb(2, vi);
    }

    m.f0.resize(m.nf); m.f1.resize(m.nf); m.f2.resize(m.nf);
    for (int fi = 0; fi < m.nf; fi++) {
        m.f0[fi] = it(0, fi) - 1;
        m.f1[fi] = it(1, fi) - 1;
        m.f2[fi] = it(2, fi) - 1;
    }

    m.nx.assign(m.nv, 0); m.ny.assign(m.nv, 0); m.nz.assign(m.nv, 0);
    m.dx.assign(m.nv, 0); m.dy.assign(m.nv, 0); m.dz.assign(m.nv, 0);
    m.odx.assign(m.nv, 0); m.ody.assign(m.nv, 0); m.odz.assign(m.nv, 0);
}

// ---------------------------------------------------------------------------
// The surface-relaxation algorithms implemented here (mris_inflate,
// mris_smooth, mris_sphere) are only ever intended to run on closed
// (watertight), manifold, genus-0 input: they are *not* robust to surfaces
// with boundaries, handles, or non-manifold edges. Vertices on a
// boundary/defect have incomplete, asymmetric 1-/2-ring neighborhoods, which
// introduces a persistent directional bias into the distance/spring forces;
// under momentum integration this compounds into a strongly distorted
// ("sausage"-shaped) result instead of converging towards a smooth surface.
//
// Isosurfaces produced by marching-cubes-style routines (e.g. vcg_isosurface)
// commonly contain small extraction cracks that violate this precondition, so
// we fail fast with an actionable message rather than silently producing
// garbage.
// ---------------------------------------------------------------------------
void mrisCheckClosedManifold(const MrisMesh &m, const char *caller)
{
    std::unordered_map<int64_t, int> edge_count;
    edge_count.reserve((size_t)m.nf * 3);
    auto add_edge = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        int64_t key = (int64_t)a * (int64_t)m.nv + (int64_t)b;
        edge_count[key]++;
    };
    for (int fi = 0; fi < m.nf; fi++) {
        add_edge(m.f0[fi], m.f1[fi]);
        add_edge(m.f1[fi], m.f2[fi]);
        add_edge(m.f2[fi], m.f0[fi]);
    }
    int boundary = 0, nonmanifold = 0;
    for (const auto &kv : edge_count) {
        if (kv.second == 1) boundary++;
        else if (kv.second > 2) nonmanifold++;
    }
    if (boundary > 0 || nonmanifold > 0) {
        Rcpp::stop(
            "%s: input mesh is not a closed, manifold surface "
            "(boundary edges: %d, non-manifold edges: %d). This algorithm "
            "requires a topologically-correct, closed (watertight), genus-0 "
            "surface. Surfaces extracted by isosurface/marching-cubes "
            "routines (e.g. vcg_isosurface) often contain small extraction "
            "cracks that violate this. "
            "Please repair the mesh (fill holes, remove non-manifold edges "
            "(e.g. via vcg_fix_defects) before calling %s().",
            caller, boundary, nonmanifold, caller
        );
    }
}

// ---------------------------------------------------------------------------
// Build 1-ring and 2-ring adjacency from the face list.
// Uses a mark[] array with per-vertex "generation" indices to avoid
// clearing the array for each vertex (O(nv * deg) total).
// ---------------------------------------------------------------------------
void mrisBuildAdjacency(MrisMesh &m, bool two_ring)
{
    // 1-ring: collect all edge-sharing neighbors, then deduplicate.
    std::vector<std::vector<int>> adj(m.nv);
    for (int fi = 0; fi < m.nf; fi++) {
        int a = m.f0[fi], b = m.f1[fi], c = m.f2[fi];
        adj[a].push_back(b); adj[a].push_back(c);
        adj[b].push_back(a); adj[b].push_back(c);
        adj[c].push_back(a); adj[c].push_back(b);
    }
    m.nbrs1.resize(m.nv);
    for (int vi = 0; vi < m.nv; vi++) {
        auto &a = adj[vi];
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
        m.nbrs1[vi] = std::move(a);
    }

    if (!two_ring) return;

    // 2-ring: union of 1-ring and neighbors-of-1-ring, excluding self.
    // mark[vj] == vi  means vj has already been added for vertex vi.
    m.nbrs2.resize(m.nv);
    std::vector<int> mark(m.nv, -1);
    for (int vi = 0; vi < m.nv; vi++) {
        mark[vi] = vi;                    // exclude self
        auto &nb2 = m.nbrs2[vi];
        for (int nj : m.nbrs1[vi]) {     // 1-ring
            if (mark[nj] != vi) { mark[nj] = vi; nb2.push_back(nj); }
        }
        for (int nj : m.nbrs1[vi]) {     // 2-ring via each 1-ring neighbor
            for (int nk : m.nbrs1[nj]) {
                if (mark[nk] != vi) { mark[nk] = vi; nb2.push_back(nk); }
            }
        }
    }

    // Average 2-ring neighborhood size, used to normalize gradient terms.
    double tot = 0.0;
    for (int vi = 0; vi < m.nv; vi++) tot += m.nbrs2[vi].size();
    m.avg_nbrs = (m.nv > 0) ? (float)(tot / m.nv) : 6.0f;

    m.dist_orig.resize(m.nv);
    m.dist_curr.resize(m.nv);
}

// ---------------------------------------------------------------------------
void mrisComputeNormals(MrisMesh &m)
{
    std::fill(m.nx.begin(), m.nx.end(), 0.0f);
    std::fill(m.ny.begin(), m.ny.end(), 0.0f);
    std::fill(m.nz.begin(), m.nz.end(), 0.0f);

    for (int fi = 0; fi < m.nf; fi++) {
        int i0 = m.f0[fi], i1 = m.f1[fi], i2 = m.f2[fi];
        float ax = m.x[i1]-m.x[i0], ay = m.y[i1]-m.y[i0], az = m.z[i1]-m.z[i0];
        float bx = m.x[i2]-m.x[i0], by = m.y[i2]-m.y[i0], bz = m.z[i2]-m.z[i0];
        // cross product = face normal, magnitude = 2 * face area
        float cx = ay*bz - az*by, cy = az*bx - ax*bz, cz = ax*by - ay*bx;
        m.nx[i0] += cx; m.ny[i0] += cy; m.nz[i0] += cz;
        m.nx[i1] += cx; m.ny[i1] += cy; m.nz[i1] += cz;
        m.nx[i2] += cx; m.ny[i2] += cy; m.nz[i2] += cz;
    }
    for (int vi = 0; vi < m.nv; vi++) {
        float l = std::sqrt(m.nx[vi]*m.nx[vi] + m.ny[vi]*m.ny[vi] + m.nz[vi]*m.nz[vi]);
        if (l > 1e-10f) { m.nx[vi] /= l; m.ny[vi] /= l; m.nz[vi] /= l; }
    }
}

// ---------------------------------------------------------------------------
void mrisComputeFaceGeometry(const MrisMesh &m, std::vector<FaceGeometry> &fg)
{
    fg.resize(m.nf);
    for (int fi = 0; fi < m.nf; fi++) {
        int i0 = m.f0[fi], i1 = m.f1[fi], i2 = m.f2[fi];
        float ax = m.x[i1]-m.x[i0], ay = m.y[i1]-m.y[i0], az = m.z[i1]-m.z[i0];
        float bx = m.x[i2]-m.x[i0], by = m.y[i2]-m.y[i0], bz = m.z[i2]-m.z[i0];
        float cx = ay*bz - az*by, cy = az*bx - ax*bz, cz = ax*by - ay*bx;
        float mag = std::sqrt(cx*cx + cy*cy + cz*cz);
        float fa  = 0.5f * mag;

        FaceGeometry &f = fg[fi];
        if (mag > 1e-10f) { f.nx = cx/mag; f.ny = cy/mag; f.nz = cz/mag; }
        else              { f.nx = 0.0f;   f.ny = 0.0f;   f.nz = 0.0f;  }

        // If the raw cross-product normal points toward the surface's center
        // (assumed at the origin; see mrisCenter), the face is inverted:
        // negate its area and flip its cached normal so it points outward.
        // This corrected, outward-pointing normal is what the area term's
        // gradient subsequently uses.
        float dot = (m.x[i0]+m.x[i1]+m.x[i2]) * f.nx
                  + (m.y[i0]+m.y[i1]+m.y[i2]) * f.ny
                  + (m.z[i0]+m.z[i1]+m.z[i2]) * f.nz;
        if (dot < 0.0f) {
            f.area = -fa;
            f.nx = -f.nx; f.ny = -f.ny; f.nz = -f.nz;
        } else {
            f.area = fa;
        }
    }
}

// ---------------------------------------------------------------------------
void mrisComputeArea(MrisMesh &m)
{
    std::vector<FaceGeometry> fg;
    mrisComputeFaceGeometry(m, fg);

    m.total_area = 0.0f;
    m.neg_area   = 0.0f;
    for (int fi = 0; fi < m.nf; fi++) {
        if (fg[fi].area >= 0.0f) m.total_area += fg[fi].area;
        else                     m.neg_area   += -fg[fi].area;
    }
}

// ---------------------------------------------------------------------------
void mrisComputeCurrentDistances(MrisMesh &m)
{
    for (int vi = 0; vi < m.nv; vi++) {
        const auto &nb = m.nbrs2[vi];
        auto &dc = m.dist_curr[vi];
        int nn = (int)nb.size();
        dc.resize(nn);
        for (int ni = 0; ni < nn; ni++) {
            int nj = nb[ni];
            float ddx = m.x[nj]-m.x[vi], ddy = m.y[nj]-m.y[vi], ddz = m.z[nj]-m.z[vi];
            dc[ni] = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
        }
    }
}

// ---------------------------------------------------------------------------
void mrisStoreOriginalDistances(MrisMesh &m)
{
    for (int vi = 0; vi < m.nv; vi++) {
        const auto &nb = m.nbrs2[vi];
        int nn = (int)nb.size();
        m.dist_orig[vi].resize(nn);
        m.dist_curr[vi].resize(nn);
        for (int ni = 0; ni < nn; ni++) {
            int nj = nb[ni];
            float ddx = m.x[nj]-m.x[vi], ddy = m.y[nj]-m.y[vi], ddz = m.z[nj]-m.z[vi];
            float d = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
            m.dist_orig[vi][ni] = d;
            m.dist_curr[vi][ni] = d;   // identical at start
        }
    }
}

// ---------------------------------------------------------------------------
void mrisScaleToArea(MrisMesh &m, float target_area)
{
    if (m.total_area <= 0.0f) return;
    float s = std::sqrt(target_area / m.total_area);
    for (int vi = 0; vi < m.nv; vi++) {
        m.x[vi] *= s; m.y[vi] *= s; m.z[vi] *= s;
    }
    // area scales as s^2, so update without recomputing
    m.neg_area   *= (s * s);
    m.total_area  = target_area;
}

// ---------------------------------------------------------------------------
void mrisCenter(MrisMesh &m)
{
    double cx = 0, cy = 0, cz = 0;
    for (int vi = 0; vi < m.nv; vi++) { cx += m.x[vi]; cy += m.y[vi]; cz += m.z[vi]; }
    float fcx = (float)(cx / m.nv), fcy = (float)(cy / m.nv), fcz = (float)(cz / m.nv);
    for (int vi = 0; vi < m.nv; vi++) { m.x[vi] -= fcx; m.y[vi] -= fcy; m.z[vi] -= fcz; }
}

// ---------------------------------------------------------------------------
// Per-vertex radial projection onto the sphere:
//   d = 1 - r / |p|;  p <- p - d * p   (i.e. p <- p * r / |p|)
// ---------------------------------------------------------------------------
void mrisProjectOntoSphere(MrisMesh &m, float r)
{
    for (int vi = 0; vi < m.nv; vi++) {
        float x = m.x[vi], y = m.y[vi], z = m.z[vi];
        float dist = std::sqrt(x*x + y*y + z*z);
        float d = (dist > 1e-10f) ? (1.0f - r / dist) : 0.0f;
        m.x[vi] = x - d * x;
        m.y[vi] = y - d * y;
        m.z[vi] = z - d * z;
    }
}

// ---------------------------------------------------------------------------
// Shared "include-self" 1-ring averaging step: `n_averages` passes applied
// to the force gradient (mris_inflate / mris_sphere) or to vertex positions
// (mris_smooth).
// ---------------------------------------------------------------------------
void mrisAverageField(const MrisMesh &m, std::vector<float> &fx,
                      std::vector<float> &fy, std::vector<float> &fz,
                      int n_averages)
{
    if (n_averages <= 0) return;

    std::vector<float> tx(m.nv), ty(m.nv), tz(m.nv);
    for (int pass = 0; pass < n_averages; pass++) {
        for (int vi = 0; vi < m.nv; vi++) {
            const auto &nb = m.nbrs1[vi];
            int nn = (int)nb.size();
            float sx = fx[vi], sy = fy[vi], sz = fz[vi];
            for (int nj : nb) { sx += fx[nj]; sy += fy[nj]; sz += fz[nj]; }
            float inv = 1.0f / (nn + 1);
            tx[vi] = sx * inv; ty[vi] = sy * inv; tz[vi] = sz * inv;
        }
        fx = tx; fy = ty; fz = tz;
    }
}

// ---------------------------------------------------------------------------
// Apply momentum to the gradient accumulator, clamp the per-vertex
// displacement, and move the vertices (no further gradient averaging applied
// at this step).
// ---------------------------------------------------------------------------
void mrisMomentumStep(MrisMesh &m, float momentum, float dt, float max_disp_mm)
{
    for (int vi = 0; vi < m.nv; vi++) {
        m.odx[vi] = dt * m.dx[vi] + momentum * m.odx[vi];
        m.ody[vi] = dt * m.dy[vi] + momentum * m.ody[vi];
        m.odz[vi] = dt * m.dz[vi] + momentum * m.odz[vi];

        float mag = std::sqrt(m.odx[vi]*m.odx[vi] + m.ody[vi]*m.ody[vi] + m.odz[vi]*m.odz[vi]);
        if (mag > max_disp_mm) {
            float s = max_disp_mm / mag;
            m.odx[vi] *= s; m.ody[vi] *= s; m.odz[vi] *= s;
        }
        m.x[vi] += m.odx[vi]; m.y[vi] += m.ody[vi]; m.z[vi] += m.odz[vi];
    }
}

// ---------------------------------------------------------------------------
// RMS of each vertex's height above the tangent plane of its 1-ring-neighbor
// centroid.
// ---------------------------------------------------------------------------
float mrisRmsTPHeight(const MrisMesh &m)
{
    double sum = 0.0;
    int cnt = 0;
    for (int vi = 0; vi < m.nv; vi++) {
        const auto &nb = m.nbrs1[vi];
        if (nb.empty()) continue;
        float cx = 0, cy = 0, cz = 0;
        for (int nj : nb) { cx += m.x[nj]; cy += m.y[nj]; cz += m.z[nj]; }
        float inv = 1.0f / (float)nb.size();
        cx *= inv; cy *= inv; cz *= inv;
        float nc = (m.x[vi]-cx)*m.nx[vi] + (m.y[vi]-cy)*m.ny[vi] + (m.z[vi]-cz)*m.nz[vi];
        sum += (double)nc * nc;
        cnt++;
    }
    return (cnt > 0) ? (float)std::sqrt(sum / cnt) : 0.0f;
}

// ---------------------------------------------------------------------------
Rcpp::List mrisPackMeshOutput(const MrisMesh &m)
{
    Rcpp::NumericMatrix out_vb(3, m.nv);
    Rcpp::NumericMatrix out_normals(3, m.nv);

    for (int vi = 0; vi < m.nv; vi++) {
        out_vb(0, vi) = m.x[vi]; out_vb(1, vi) = m.y[vi]; out_vb(2, vi) = m.z[vi];
        out_normals(0, vi) = m.nx[vi]; out_normals(1, vi) = m.ny[vi]; out_normals(2, vi) = m.nz[vi];
    }

    return Rcpp::List::create(
        Rcpp::Named("vb")      = out_vb,
        Rcpp::Named("normals") = out_normals
    );
}

} // namespace ravetools
