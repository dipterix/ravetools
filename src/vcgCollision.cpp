// ---------------------------------------------------------------------------
// Exact proximity / collision detection between two geometries.
//
// One side (`x`, the ROI) is indexed once into a uniform grid; the other side
// (`y`) is streamed against it. Either side may be a point cloud, a set of
// polylines, or a triangular mesh, giving nine combinations that all share a
// single code path: build the query primitive's bounding box, inflate it by
// `radius`, collect the grid cells it touches, then measure the survivors
// exactly.
//
// The result is exact in both directions: the inflated box is a conservative
// superset of everything within `radius` (no false negatives) and every
// candidate is then measured with an exact primitive distance (no false
// positives).
//
// Note on the distance primitives: vcglib's SegmentSegmentDistance
// (vcg/space/distance3.h) clamps the two infinite-line closest points onto
// their own segments *independently*, which overestimates whenever the
// clamped point on one segment is not mutually closest to the clamped point
// on the other. Measured on random skew pairs it overestimates ~44% of the
// time, by up to 3.9x. TriangleSegmentDistance and TriangleTriangleDistance
// both call it and inherit the error. Those three are therefore reimplemented
// here (segment_segment_distance uses the standard clamped algorithm, and the
// triangle cases compose it); the remaining primitives are taken from vcglib,
// which computes them exactly.
// ---------------------------------------------------------------------------

#include "vcgCommon.h"

#include <vcg/space/distance3.h>
#include <vcg/space/intersection3.h>
#include <vcg/complex/algorithms/inside.h>

#include <atomic>
#include <cmath>
#include <limits>
#include <vector>

#include "TinyParallel.h"

namespace ravetools {
namespace collision {

typedef float                Scalar;
typedef vcg::Point3<Scalar>  Point3;
typedef vcg::Box3<Scalar>    Box3;
typedef vcg::Segment3<Scalar> Segment3;
typedef vcg::Triangle3<Scalar> Triangle3;

static const Scalar kEps = std::numeric_limits<Scalar>::epsilon();

// Modes, shared with the R wrapper.
enum Mode { MODE_POINTS = 0, MODE_SEGMENTS = 1, MODE_MESH = 2 };

// How exhaustively the target is scanned, shared with the R wrapper. Results
// are always reported per *unit* of the target -- a point, a polyline, or a
// face -- and the level decides which element of a hit unit is reported and how
// early the scan gives up:
//
//   ELEMENT  every element of every unit; reports the closest one
//   UNIT     stops at the first hit inside each unit; reports that one
//   WHOLE    stops at the first hit anywhere; reports that single witness
//
// A point and a face each *are* their own unit, so ELEMENT and UNIT differ only
// for polylines.
enum Level { LEVEL_ELEMENT = 0, LEVEL_UNIT = 1, LEVEL_WHOLE = 2 };

// ---------------------------------------------------------------------------
// Exact primitive distances
// ---------------------------------------------------------------------------

inline Scalar point_point_distance(const Point3 &p, const Point3 &q)
{
  return (p - q).Norm();
}

inline Scalar segment_point_distance(const Segment3 &s, const Point3 &p)
{
  Point3 closest;
  Scalar sqr_dist;
  vcg::SegmentPointSquaredDistance(s, p, closest, sqr_dist);
  if (sqr_dist <= 0) { return 0; }
  return std::sqrt(sqr_dist);
}

// Closest distance between two segments, clamping both parameters jointly.
// Handles degenerate (zero length) segments, which the vcglib version does
// not: it normalizes the segment directions and yields NaN for those.
inline Scalar segment_segment_distance(const Point3 &p0, const Point3 &q0,
                                       const Point3 &p1, const Point3 &q1)
{
  const Point3 d0 = q0 - p0;
  const Point3 d1 = q1 - p1;
  const Point3 r  = p0 - p1;

  const Scalar a = d0.dot(d0);   // squared length of segment 0
  const Scalar e = d1.dot(d1);   // squared length of segment 1
  const Scalar f = d1.dot(r);

  Scalar s, t;

  if (a <= kEps && e <= kEps) {
    // Both segments degenerate to points.
    return r.Norm();
  }
  if (a <= kEps) {
    // Segment 0 is a point.
    s = 0;
    t = f / e;
    t = (t < 0) ? 0 : ((t > 1) ? 1 : t);
  } else {
    const Scalar c = d0.dot(r);
    if (e <= kEps) {
      // Segment 1 is a point.
      t = 0;
      s = -c / a;
      s = (s < 0) ? 0 : ((s > 1) ? 1 : s);
    } else {
      const Scalar b     = d0.dot(d1);
      const Scalar denom = a * e - b * b;

      // denom == 0 means the segments are parallel; any s works, so pick 0
      // and let the t clamp below resolve it.
      if (denom > kEps) {
        s = (b * f - c * e) / denom;
        s = (s < 0) ? 0 : ((s > 1) ? 1 : s);
      } else {
        s = 0;
      }

      t = (b * s + f) / e;

      // Clamping t invalidates s, so recompute s for the clamped t. This
      // step is what vcglib's version omits.
      if (t < 0) {
        t = 0;
        s = -c / a;
        s = (s < 0) ? 0 : ((s > 1) ? 1 : s);
      } else if (t > 1) {
        t = 1;
        s = (b - c) / a;
        s = (s < 0) ? 0 : ((s > 1) ? 1 : s);
      }
    }
  }

  const Point3 c0 = p0 + d0 * s;
  const Point3 c1 = p1 + d1 * t;
  return (c0 - c1).Norm();
}

inline Scalar segment_segment_distance(const Segment3 &s0, const Segment3 &s1)
{
  return segment_segment_distance(s0.P0(), s0.P1(), s1.P0(), s1.P1());
}

inline Scalar triangle_point_distance(const Point3 &v0, const Point3 &v1,
                                      const Point3 &v2, const Point3 &p)
{
  Scalar dist;
  Point3 closest;
  vcg::TrianglePointDistance(Triangle3(v0, v1, v2), p, dist, closest);
  return dist;
}

// Minimum distance between a triangle and a segment.
//
// When they do not intersect, the closest pair is either (segment endpoint,
// triangle) or (segment interior, triangle edge): an interior-to-interior
// optimum requires the segment to be parallel to the triangle plane, in which
// case the endpoints attain the same value. So the five terms below are
// exhaustive.
inline Scalar triangle_segment_distance(const Point3 &v0, const Point3 &v1,
                                        const Point3 &v2, const Segment3 &s)
{
  Scalar a, b;
  if (vcg::IntersectionSegmentTriangle(s, v0, v1, v2, a, b)) { return 0; }

  Scalar best = triangle_point_distance(v0, v1, v2, s.P0());

  Scalar d = triangle_point_distance(v0, v1, v2, s.P1());
  if (d < best) { best = d; }

  d = segment_segment_distance(s.P0(), s.P1(), v0, v1);
  if (d < best) { best = d; }
  d = segment_segment_distance(s.P0(), s.P1(), v1, v2);
  if (d < best) { best = d; }
  d = segment_segment_distance(s.P0(), s.P1(), v2, v0);
  if (d < best) { best = d; }

  return best;
}

// Minimum distance between two triangles: the closest pair always involves at
// least one edge, so testing all six edges against the opposite triangle is
// exhaustive. Intersecting triangles always have an edge crossing the other,
// which triangle_segment_distance reports as 0.
inline Scalar triangle_triangle_distance(const Point3 &a0, const Point3 &a1,
                                         const Point3 &a2, const Point3 &b0,
                                         const Point3 &b1, const Point3 &b2)
{
  Scalar best = triangle_segment_distance(b0, b1, b2, Segment3(a0, a1));
  if (best <= 0) { return 0; }

  Scalar d = triangle_segment_distance(b0, b1, b2, Segment3(a1, a2));
  if (d < best) { best = d; }
  if (best <= 0) { return 0; }
  d = triangle_segment_distance(b0, b1, b2, Segment3(a2, a0));
  if (d < best) { best = d; }
  if (best <= 0) { return 0; }

  d = triangle_segment_distance(a0, a1, a2, Segment3(b0, b1));
  if (d < best) { best = d; }
  if (best <= 0) { return 0; }
  d = triangle_segment_distance(a0, a1, a2, Segment3(b1, b2));
  if (d < best) { best = d; }
  if (best <= 0) { return 0; }
  d = triangle_segment_distance(a0, a1, a2, Segment3(b2, b0));
  if (d < best) { best = d; }

  return best;
}

// ---------------------------------------------------------------------------
// Grid elements
//
// GridStaticPtr<OBJTYPE, FLT> asks only three things of OBJTYPE: a ScalarType
// typedef, GetBBox(Box3<ScalarType>&), and IsD(). Mesh ROIs use MyMesh's own
// face type so the very same grid can also drive vcg::tri::Inside::Is_Inside;
// point and polyline ROIs use the two small structs below.
// ---------------------------------------------------------------------------

struct PointElem {
  typedef Scalar ScalarType;
  Point3 p;
  int    idx;   // 0-based row of the source matrix

  inline void GetBBox(Box3 &bb) const { bb.min = p; bb.max = p; }
  inline bool IsD() const { return false; }
};

struct SegElem {
  typedef Scalar ScalarType;
  Point3 a, b;
  int    idx;   // 0-based row where this segment starts

  inline void GetBBox(Box3 &bb) const
  {
    bb.SetNull();
    bb.Add(a);
    bb.Add(b);
  }
  inline bool IsD() const { return false; }
};

// vcg's FaceTmark mutates per-face IMark, which races across threads. Queries
// here are read-only apart from that marker, so a no-op marker makes the whole
// query thread-safe. The cost is that an element straddling several cells may
// be measured more than once, which changes nothing but a little work.
struct NoopMark {
  inline void UnMarkAll() const {}
  template <typename T> inline bool IsMarked(T *) const { return false; }
  template <typename T> inline void Mark(T *) const {}
  inline void SetMesh(void *) const {}
};

// Source row of an ROI element, for the reported `index`.
inline int elem_row(const PointElem &e, MyMesh *) { return e.idx; }
inline int elem_row(const SegElem &e, MyMesh *) { return e.idx; }
inline int elem_row(const MyFace &f, MyMesh *m)
{
  return static_cast<int>(vcg::tri::Index(*m, &f));
}

// Bounding box of an ROI element.
inline void elem_bbox(const PointElem &e, Box3 &bb) { e.GetBBox(bb); }
inline void elem_bbox(const SegElem &e, Box3 &bb) { e.GetBBox(bb); }
inline void elem_bbox(const MyFace &f, Box3 &bb) { f.GetBBox(bb); }

// Distance from an ROI element to each of the three query primitives.
inline Scalar elem_distance(const PointElem &e, const Point3 &q)
{
  return point_point_distance(e.p, q);
}
inline Scalar elem_distance(const PointElem &e, const Segment3 &s)
{
  return segment_point_distance(s, e.p);
}
inline Scalar elem_distance(const PointElem &e, const Triangle3 &t)
{
  return triangle_point_distance(t.cP(0), t.cP(1), t.cP(2), e.p);
}

inline Scalar elem_distance(const SegElem &e, const Point3 &q)
{
  return segment_point_distance(Segment3(e.a, e.b), q);
}
inline Scalar elem_distance(const SegElem &e, const Segment3 &s)
{
  return segment_segment_distance(e.a, e.b, s.P0(), s.P1());
}
inline Scalar elem_distance(const SegElem &e, const Triangle3 &t)
{
  return triangle_segment_distance(t.cP(0), t.cP(1), t.cP(2), Segment3(e.a, e.b));
}

inline Scalar elem_distance(const MyFace &f, const Point3 &q)
{
  return triangle_point_distance(f.cP(0), f.cP(1), f.cP(2), q);
}
inline Scalar elem_distance(const MyFace &f, const Segment3 &s)
{
  return triangle_segment_distance(f.cP(0), f.cP(1), f.cP(2), s);
}
inline Scalar elem_distance(const MyFace &f, const Triangle3 &t)
{
  return triangle_triangle_distance(f.cP(0), f.cP(1), f.cP(2),
                                    t.cP(0), t.cP(1), t.cP(2));
}

// ---------------------------------------------------------------------------
// The ROI index
// ---------------------------------------------------------------------------

// Result of one proximity query.
struct Hit {
  bool   found;
  Scalar distance;
  int    index;   // 0-based ROI row, or -1
  int    sub;     // 1-based element within the target unit

  Hit() : found(false), distance(0), index(-1), sub(1) {}
};

// Which element of the query primitive to name in the result. A point and a
// segment are atomic, so the answer is 1 (for a segment the driver overwrites
// it with the ordinal inside the polyline, which only the driver knows). A
// triangle reports its vertex nearest the ROI element -- a locator only: the
// closest point on a triangle is often on an edge or in the interior, and the
// reported distance stays the exact triangle distance.
template <class ElemType>
inline int sub_index(const Point3 &, const ElemType &) { return 1; }

template <class ElemType>
inline int sub_index(const Segment3 &, const ElemType &) { return 1; }

template <class ElemType>
inline int sub_index(const Triangle3 &t, const ElemType &e)
{
  int    best = 0;
  Scalar bd   = elem_distance(e, t.cP(0));

  Scalar d = elem_distance(e, t.cP(1));
  if (d < bd) { bd = d; best = 1; }
  d = elem_distance(e, t.cP(2));
  if (d < bd) { best = 2; }

  return best + 1;
}

// Scratch buffer reused across queries by one thread.
template <class ElemType>
struct QueryScratch {
  std::vector<ElemType *> candidates;
};

template <class ElemType>
class CollisionIndex {
public:
  typedef vcg::GridStaticPtr<ElemType, Scalar> GridType;
  typedef QueryScratch<ElemType>               Scratch;

  CollisionIndex() : radius_(0), chunk_(1), pad_(0), mesh_(NULL), empty_(true) {}

  // `first`/`last` must stay alive and pointer-stable for the index's lifetime.
  template <class ObjIter>
  void build(ObjIter first, ObjIter last, std::size_t n, Scalar radius,
             MyMesh *mesh)
  {
    radius_ = radius;
    mesh_   = mesh;
    empty_  = (n == 0);
    if (empty_) { return; }

    // Bounding box of the whole ROI, used both for the O(1) reject and to
    // pick a cell size.
    Box3 bb;
    bb.SetNull();
    for (ObjIter it = first; it != last; ++it) {
      Box3 eb;
      elem_bbox(*it, eb);
      bb.Add(eb);
    }

    // Cell size. The default is "about as many cells as elements"; enlarge it
    // to `radius` when the radius is coarser, so a query box spans O(1) cells
    // rather than O((radius / cell)^3). Then floor it so the grid cannot blow
    // up in memory on a tiny radius.
    const Scalar diag = bb.Diag();
    Scalar cell = radius;
    if (diag > 0) {
      const Scalar automatic =
        static_cast<Scalar>(diag / std::max(1.0, std::cbrt(static_cast<double>(n))));
      if (automatic > cell) { cell = automatic; }
      const Scalar floor_cell = static_cast<Scalar>(diag / 200.0);  // <= 200^3 cells
      if (cell < floor_cell) { cell = floor_cell; }
    }
    if (!(cell > 0)) { cell = 1; }

    grid_.SetWithRadius(first, last, cell);

    roi_box_ = grid_.bbox;

    // vcg's Box3::Collide uses strict inequalities, so a box with zero extent
    // on any axis never collides -- which is exactly what a point query, an
    // axis-aligned segment, or a radius of 0 produces. Widening the *search*
    // box by a hair restores those cases and is always safe: it can only admit
    // extra candidates, which the exact distance test then rejects.
    pad_ = std::max<Scalar>(diag * static_cast<Scalar>(1e-6), static_cast<Scalar>(1e-6));

    // Query boxes are kept to O(1) cells by splitting long segments into
    // chunks of roughly one cell.
    chunk_ = std::max(std::max(grid_.voxel[0], grid_.voxel[1]), grid_.voxel[2]);
    if (!(chunk_ > 0)) { chunk_ = cell; }
  }

  bool empty() const { return empty_; }
  Scalar radius() const { return radius_; }
  MyMesh *mesh() const { return mesh_; }
  GridType &grid() { return grid_; }

  Hit test_point(const Point3 &p, Scratch &scratch)
  {
    Box3 box;
    box.min = p;
    box.max = p;
    return query(box, p, scratch);
  }

  Hit test_segment(const Point3 &a, const Point3 &b, Scratch &scratch)
  {
    Hit best;

    const Point3 d   = b - a;
    const Scalar len = d.Norm();

    // Split long segments so each query box stays small. Exact: the minimum
    // over the chunks is the minimum over the whole segment.
    int nchunk = 1;
    if (chunk_ > 0 && len > chunk_) {
      double raw = std::ceil(static_cast<double>(len / chunk_));
      if (raw > 4096.0) { raw = 4096.0; }
      nchunk = static_cast<int>(raw);
    }

    for (int k = 0; k < nchunk; k++) {
      const Scalar t0 = static_cast<Scalar>(k) / static_cast<Scalar>(nchunk);
      const Scalar t1 = static_cast<Scalar>(k + 1) / static_cast<Scalar>(nchunk);
      const Point3 p0 = a + d * t0;
      const Point3 p1 = a + d * t1;

      Box3 box;
      box.SetNull();
      box.Add(p0);
      box.Add(p1);

      const Hit h = query(box, Segment3(p0, p1), scratch);
      if (h.found && (!best.found || h.distance < best.distance)) { best = h; }
      if (best.found && best.distance <= 0) { break; }
    }
    return best;
  }

  Hit test_triangle(const Point3 &v0, const Point3 &v1, const Point3 &v2,
                    Scratch &scratch)
  {
    Box3 box;
    box.SetNull();
    box.Add(v0);
    box.Add(v1);
    box.Add(v2);
    return query(box, Triangle3(v0, v1, v2), scratch);
  }

  // Whole-polyline / whole-mesh reject: does anything in this box come within
  // `radius` of the ROI at all?
  bool box_may_collide(const Box3 &box) const
  {
    if (empty_) { return false; }
    Box3 inflated = box;
    inflated.Offset(radius_ + pad_);
    return inflated.Collide(roi_box_);
  }

private:
  template <class Primitive>
  Hit query(Box3 box, const Primitive &prim, Scratch &scratch)
  {
    Hit result;
    if (empty_) { return result; }

    box.Offset(radius_ + pad_);
    if (!box.Collide(roi_box_)) { return result; }

    NoopMark marker;
    grid_.GetInBox(marker, box, scratch.candidates);

    const std::size_t n = scratch.candidates.size();
    for (std::size_t i = 0; i < n; i++) {
      ElemType *elem = scratch.candidates[i];
      const Scalar d = elem_distance(*elem, prim);
      if (d <= radius_ && (!result.found || d < result.distance)) {
        result.found    = true;
        result.distance = d;
        result.index    = elem_row(*elem, mesh_);
        result.sub      = sub_index(prim, *elem);
      }
    }
    return result;
  }

  GridType grid_;
  Box3     roi_box_;
  Scalar   radius_;
  Scalar   chunk_;
  Scalar   pad_;
  MyMesh  *mesh_;
  bool     empty_;
};

// ---------------------------------------------------------------------------
// Target side
// ---------------------------------------------------------------------------

// A run of consecutive, non-separator rows. May be empty.
struct Group {
  int begin;
  int end;   // exclusive
};

// Splits a 3 x n coordinate buffer into the chains delimited by NaN columns.
//
// Every slot between separators is a chain, empty ones included, so `A NA NA B`
// is three chains and the middle one is empty. The two ends are not symmetric:
// a separator in the first row opens an empty leading chain, but a separator in
// the last row terminates the final chain rather than opening a trailing empty
// one. That is what makes the idiom
//
//   do.call("rbind", lapply(tracts, function(line) rbind(line, NA)))
//
// yield exactly one chain per tract -- an empty tract contributes only its own
// separator, and the separator closing the last tract adds nothing.
inline std::vector<Group> find_groups(const std::vector<Scalar> &coords, int n)
{
  std::vector<Group> groups;
  if (n == 0) { return groups; }

  int start = 0;
  for (int i = 0; i < n; i++) {
    if (coords[i * 3] == coords[i * 3]) { continue; }   // not a separator
    Group g = { start, i };
    groups.push_back(g);
    start = i + 1;
  }
  // Skipped when the last row is itself a separator, which closes the chain
  // before it instead of starting an empty one.
  if (start < n) {
    Group g = { start, n };
    groups.push_back(g);
  }
  return groups;
}

inline Point3 coord_at(const std::vector<Scalar> &coords, int i)
{
  return Point3(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]);
}

// One slot per unit of the target, in plain buffers so worker threads never
// touch the R API. Each group of the target owns its own slots, so workers
// never write to the same slot and no synchronization is needed.
struct UnitOutput {
  std::vector<int>    hit;        // 1 = hit, 0 = miss, NA_INTEGER = nothing to test
  std::vector<int>    index;      // 1-based element within the unit
  std::vector<double> distance;
  std::vector<int>    x_index;    // 1-based ROI element, NA_INTEGER when unknown

  // Slots start at NA, and a unit with nothing to test is simply never written:
  // a separator row, an empty or single-vertex chain, a deleted face. That NA is
  // part of the contract and reaches R untouched -- it distinguishes "no
  // geometry here" from the 0 that `record` writes for "tested and missed".
  void resize(std::size_t n)
  {
    hit.assign(n, NA_INTEGER);
    index.assign(n, NA_INTEGER);
    distance.assign(n, NA_REAL);
    x_index.assign(n, NA_INTEGER);
  }

  void record(std::size_t i, const Hit &h)
  {
    if (h.found) {
      hit[i]      = 1;
      index[i]    = h.sub;
      distance[i] = static_cast<double>(h.distance);
      x_index[i]  = (h.index >= 0) ? (h.index + 1) : NA_INTEGER;
    } else {
      hit[i] = 0;
    }
  }
};

// The single witness recorded at LEVEL_WHOLE. `flag` is the short-circuit
// signal every worker polls; the first thread to flip it is the only one that
// writes the payload, so no mutex is needed. Which unit wins is therefore
// timing-dependent when several threads are running -- `flag` itself is not.
struct WholeResult {
  std::atomic<bool> flag;
  int               unit;       // 0-based unit ordinal
  int               sub;
  double            distance;
  int               x_index;    // 1-based, or NA_INTEGER

  WholeResult()
    : flag(false), unit(-1), sub(1), distance(NA_REAL), x_index(NA_INTEGER) {}

  bool hit() const { return flag.load(); }

  void claim(int unit_, const Hit &h)
  {
    if (flag.exchange(true)) { return; }   // someone else got there first
    unit     = unit_;
    sub      = h.sub;
    distance = static_cast<double>(h.distance);
    x_index  = (h.index >= 0) ? (h.index + 1) : NA_INTEGER;
  }
};

// Optional "is this point strictly inside the closed ROI surface" test. Only
// meaningful for a mesh ROI, and only reached when the radius test missed.
template <class ElemType>
struct InteriorTester {
  // Non-mesh ROI: never applies.
  static bool inside(CollisionIndex<ElemType> &, const Point3 &) { return false; }
};

template <>
struct InteriorTester<MyFace> {
  static bool inside(CollisionIndex<MyFace> &index, const Point3 &p)
  {
    MyMesh *m = index.mesh();
    if (m == NULL || m->fn == 0) { return false; }
    if (!m->bbox.IsIn(p)) { return false; }
    typedef CollisionIndex<MyFace>::GridType GridType;
    return vcg::tri::Inside<GridType, MyMesh>::Is_Inside(*m, index.grid(), p);
  }
};

// One unit of work: scan a group of target rows (or the target mesh's faces)
// against the index.
template <class ElemType>
struct Driver {
  typedef typename CollisionIndex<ElemType>::Scratch Scratch;

  CollisionIndex<ElemType> *index;
  const std::vector<Scalar> *coords;   // 3 x n, points / segments modes
  const std::vector<Group>  *groups;
  MyMesh                    *target_mesh;   // mesh mode
  UnitOutput                *out;           // NULL at LEVEL_WHOLE
  WholeResult               *whole;         // only used at LEVEL_WHOLE
  int                        mode_y;
  int                        level;
  bool                       include_interior;

  // True once some thread has answered the LEVEL_WHOLE question.
  inline bool answered() const
  {
    return level == LEVEL_WHOLE && whole->hit();
  }

  void run_group(std::size_t gi, Scratch &scratch) const
  {
    if (mode_y == MODE_MESH) {
      run_mesh(scratch);
      return;
    }

    const Group g = (*groups)[gi];

    // An empty chain has no vertices to box and nothing to test; its slot keeps
    // the NA it was initialized with.
    if (g.begin >= g.end) { return; }

    // Whole-polyline reject: one box test can skip an entire tract.
    Box3 gbox;
    gbox.SetNull();
    for (int i = g.begin; i < g.end; i++) { gbox.Add(coord_at(*coords, i)); }
    const bool worth_testing = index->box_may_collide(gbox);

    if (mode_y == MODE_SEGMENTS) {
      run_polyline(gi, g, worth_testing, scratch);
    } else {
      run_points(g, worth_testing, scratch);
    }
  }

  // Each row is its own unit, so its slot is the row itself. Separator rows lie
  // outside every group and keep their initial NA.
  void run_points(const Group &g, bool worth_testing, Scratch &scratch) const
  {
    for (int i = g.begin; i < g.end; i++) {
      if (answered()) { return; }

      const Point3 p = coord_at(*coords, i);

      Hit h;
      if (worth_testing) { h = index->test_point(p, scratch); }
      if (!h.found && include_interior &&
          InteriorTester<ElemType>::inside(*index, p)) {
        h.found    = true;
        h.distance = 0;
        h.index    = -1;
      }

      if (level == LEVEL_WHOLE) {
        if (h.found) { whole->claim(i, h); return; }
      } else {
        out->record(i, h);
      }
    }
  }

  // The whole polyline is one unit, occupying the slot of its group ordinal.
  // A group of a single row carries no segment at all and stays NA.
  void run_polyline(std::size_t gi, const Group &g, bool worth_testing,
                    Scratch &scratch) const
  {
    if (g.end - g.begin < 2) { return; }

    Hit best;

    for (int i = g.begin; i + 1 < g.end; i++) {
      if (answered()) { return; }

      Hit h;
      if (worth_testing) {
        h = index->test_segment(coord_at(*coords, i), coord_at(*coords, i + 1),
                                scratch);
      }
      if (!h.found && include_interior &&
          InteriorTester<ElemType>::inside(*index, coord_at(*coords, i))) {
        h.found    = true;
        h.distance = 0;
        h.index    = -1;
      }
      if (!h.found) { continue; }

      h.sub = i - g.begin + 1;   // 1-based, relative to the polyline start

      if (level == LEVEL_WHOLE) {
        whole->claim(static_cast<int>(gi), h);
        return;
      }
      if (level == LEVEL_UNIT) {
        out->record(gi, h);   // first hit wins; the rest is not measured
        return;
      }
      if (!best.found || h.distance < best.distance) { best = h; }
    }

    // Nothing to write at LEVEL_WHOLE: reaching here means this polyline missed,
    // and only the first hit anywhere is reported.
    if (level != LEVEL_WHOLE) {
      out->record(gi, best);   // records a miss when nothing was found
    }
  }

  // Each face is its own unit, so its slot is the face itself. A face cannot
  // short-circuit within itself, so ELEMENT and UNIT do the same work.
  void run_mesh(Scratch &scratch) const
  {
    const int nf = target_mesh->fn;
    for (int i = 0; i < nf; i++) {
      if (answered()) { return; }

      const MyFace &f = target_mesh->face[i];
      if (f.IsD()) { continue; }

      Hit h = index->test_triangle(f.cP(0), f.cP(1), f.cP(2), scratch);

      if (!h.found && include_interior) {
        const Point3 bary = (f.cP(0) + f.cP(1) + f.cP(2)) / 3.0f;
        if (InteriorTester<ElemType>::inside(*index, bary)) {
          h.found    = true;
          h.distance = 0;
          h.index    = -1;   // no ROI element and no nearest vertex to name
        }
      }

      if (level == LEVEL_WHOLE) {
        if (h.found) { whole->claim(i, h); return; }
      } else {
        out->record(i, h);
      }
    }
  }
};

template <class ElemType>
struct ParallelDriver : public TinyParallel::Worker {
  const Driver<ElemType> *driver;

  explicit ParallelDriver(const Driver<ElemType> *d) : driver(d) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    typename CollisionIndex<ElemType>::Scratch scratch;
    for (std::size_t gi = begin; gi < end; gi++) { driver->run_group(gi, scratch); }
  }
};

// Thread count follows the package-wide `ravetools_threads()` setting, which
// TinyParallel reads from RAVETOOLS_NUM_THREADS.
template <class ElemType>
void run_driver(const Driver<ElemType> &driver, std::size_t n_groups, bool serial)
{
  if (n_groups == 0) { return; }

  if (serial || n_groups < 2) {
    typename CollisionIndex<ElemType>::Scratch scratch;
    for (std::size_t gi = 0; gi < n_groups; gi++) {
      driver.run_group(gi, scratch);
      if ((gi & 0xFF) == 0) { Rcpp::checkUserInterrupt(); }
    }
    return;
  }

  ParallelDriver<ElemType> worker(&driver);
  TinyParallel::parallelFor(0, n_groups, worker);
}

// ---------------------------------------------------------------------------
// Input marshalling
// ---------------------------------------------------------------------------

// A 3 x n R matrix becomes a flat float buffer; non-finite columns become NaN
// separators. The R wrapper guarantees a column is either wholly finite or
// wholly NA.
inline std::vector<Scalar> read_coords(const Rcpp::NumericMatrix &m, int &n)
{
  n = static_cast<int>(m.ncol());
  std::vector<Scalar> coords(static_cast<std::size_t>(n) * 3);
  for (int i = 0; i < n; i++) {
    const double x = m(0, i), y = m(1, i), z = m(2, i);
    // NA_real_ is a NaN, so R_finite() rejects NA, NaN and +/-Inf alike.
    if (R_finite(x) && R_finite(y) && R_finite(z)) {
      coords[i * 3]     = static_cast<Scalar>(x);
      coords[i * 3 + 1] = static_cast<Scalar>(y);
      coords[i * 3 + 2] = static_cast<Scalar>(z);
    } else {
      const Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();
      coords[i * 3] = coords[i * 3 + 1] = coords[i * 3 + 2] = nan;
    }
  }
  return coords;
}

inline std::vector<PointElem> build_point_elems(const std::vector<Scalar> &coords,
                                                int n)
{
  std::vector<PointElem> elems;
  elems.reserve(n);
  for (int i = 0; i < n; i++) {
    if (!(coords[i * 3] == coords[i * 3])) { continue; }   // separator
    PointElem e;
    e.p   = coord_at(coords, i);
    e.idx = i;
    elems.push_back(e);
  }
  return elems;
}

// `groups` comes from find_groups; the caller keeps it so that the chain count
// reported in the summary and the elements indexed here share one grouping.
inline std::vector<SegElem> build_seg_elems(const std::vector<Scalar> &coords,
                                            const std::vector<Group> &groups)
{
  std::vector<SegElem> elems;
  for (std::size_t g = 0; g < groups.size(); g++) {
    for (int i = groups[g].begin; i + 1 < groups[g].end; i++) {
      SegElem e;
      e.a   = coord_at(coords, i);
      e.b   = coord_at(coords, i + 1);
      e.idx = i;
      elems.push_back(e);
    }
  }
  return elems;
}

// Runs the whole query once the ROI element type is known. The returned list is
// flat -- the R wrapper assembles `unit` / `index` / `distance` / `x_index`
// into a data frame, drops `hit_unit` at LEVEL_WHOLE, and folds the two unit
// counts into the summary.
template <class ElemType, class ObjIter>
Rcpp::List dispatch_target(ObjIter first, ObjIter last, std::size_t n_elems,
                           MyMesh *roi_mesh, double radius,
                           const std::vector<Scalar> &y_coords, int ny,
                           MyMesh *y_mesh, int mode_y, int level,
                           bool include_interior, int n_units_x)
{
  CollisionIndex<ElemType> index;
  index.build(first, last, n_elems, static_cast<Scalar>(radius), roi_mesh);

  std::vector<Group> groups;
  if (mode_y == MODE_MESH) {
    Group g = { 0, y_mesh->fn };
    groups.push_back(g);
  } else {
    groups = find_groups(y_coords, ny);
  }

  // One slot per unit of the target: a face, a polyline, or a row (separator
  // rows included, so the answers stay aligned with the input matrix).
  std::size_t n_units;
  if (mode_y == MODE_MESH) {
    n_units = static_cast<std::size_t>(y_mesh->fn);
  } else if (mode_y == MODE_SEGMENTS) {
    n_units = groups.size();
  } else {
    n_units = static_cast<std::size_t>(ny);
  }

  UnitOutput  out;
  WholeResult whole;
  if (level != LEVEL_WHOLE) { out.resize(n_units); }

  Driver<ElemType> driver;
  driver.index            = &index;
  driver.coords           = &y_coords;
  driver.groups           = &groups;
  driver.target_mesh      = y_mesh;
  driver.out              = (level == LEVEL_WHOLE) ? NULL : &out;
  driver.whole            = &whole;
  driver.mode_y           = mode_y;
  driver.level            = level;
  driver.include_interior = include_interior;

  // Is_Inside builds its own mutating FaceTmark, so the interior pass must
  // stay on one thread.
  run_driver(driver, groups.size(), include_interior);

  if (level == LEVEL_WHOLE) {
    const bool found = whole.hit();
    const int  n     = found ? 1 : 0;

    Rcpp::IntegerVector rep_unit(n), rep_index(n), rep_x(n);
    Rcpp::NumericVector rep_dist(n);
    if (found) {
      rep_unit[0]  = whole.unit + 1;
      rep_index[0] = whole.sub;
      rep_dist[0]  = whole.distance;
      rep_x[0]     = whole.x_index;
    }

    return Rcpp::List::create(
      Rcpp::Named("collide")   = found,
      Rcpp::Named("hit_unit")  = R_NilValue,
      Rcpp::Named("unit")      = rep_unit,
      Rcpp::Named("index")     = rep_index,
      Rcpp::Named("distance")  = rep_dist,
      Rcpp::Named("x_index")   = rep_x,
      Rcpp::Named("n_units_x") = n_units_x,
      Rcpp::Named("n_units_y") = static_cast<int>(n_units));
  }

  // Compact the per-unit slots into the sparse table of hit units.
  std::size_t n_hit = 0;
  for (std::size_t i = 0; i < n_units; i++) {
    if (out.hit[i] == 1) { n_hit++; }
  }

  Rcpp::LogicalVector hit_unit(n_units);
  Rcpp::IntegerVector rep_unit(n_hit), rep_index(n_hit), rep_x(n_hit);
  Rcpp::NumericVector rep_dist(n_hit);

  std::size_t k = 0;
  for (std::size_t i = 0; i < n_units; i++) {
    hit_unit[i] = (out.hit[i] == NA_INTEGER) ? NA_LOGICAL : out.hit[i];
    if (out.hit[i] != 1) { continue; }
    rep_unit[k]  = static_cast<int>(i) + 1;
    rep_index[k] = out.index[i];
    rep_dist[k]  = out.distance[i];
    rep_x[k]     = out.x_index[i];
    k++;
  }

  return Rcpp::List::create(
    Rcpp::Named("collide")   = (n_hit > 0),
    Rcpp::Named("hit_unit")  = hit_unit,
    Rcpp::Named("unit")      = rep_unit,
    Rcpp::Named("index")     = rep_index,
    Rcpp::Named("distance")  = rep_dist,
    Rcpp::Named("x_index")   = rep_x,
    Rcpp::Named("n_units_x") = n_units_x,
    Rcpp::Named("n_units_y") = static_cast<int>(n_units));
}

}   // namespace collision
}   // namespace ravetools

// ---------------------------------------------------------------------------
// R entry point
//
//   mode_x / mode_y : 0 = points, 1 = segments, 2 = mesh
//   test_level      : 0 = element, 1 = unit, 2 = whole
//   x_ / y_         : 3 x n coordinate matrices (mesh `vb` in mesh mode)
//   x_it_ / y_it_   : 3 x m 0-based face matrices, or NULL
//
// Parallelism follows ravetools_threads() (RAVETOOLS_NUM_THREADS).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP vcgDetectCollision(
    const Rcpp::NumericMatrix &x_, SEXP x_it_,
    const Rcpp::NumericMatrix &y_, SEXP y_it_,
    int mode_x, int mode_y,
    double radius,
    int test_level,
    bool include_interior)
{
  using namespace ravetools;
  using namespace ravetools::collision;

  try {
    if (!(radius >= 0)) {
      Rcpp::stop("vcgDetectCollision: `radius` must be a non-negative number.");
    }
    if (test_level < LEVEL_ELEMENT || test_level > LEVEL_WHOLE) {
      Rcpp::stop("vcgDetectCollision: `test_level` must be 0, 1, or 2.");
    }
    if (include_interior && mode_x != MODE_MESH) {
      Rcpp::stop("vcgDetectCollision: `include_interior` requires a mesh `x`.");
    }

    int ny = 0;
    std::vector<Scalar> y_coords;
    MyMesh y_mesh;

    if (mode_y == MODE_MESH) {
      if (IOMesh<MyMesh>::vcgReadR(y_mesh, Rcpp::wrap(y_), y_it_) != 0) {
        Rcpp::stop("vcgDetectCollision: `y` has no faces; cannot use mesh mode.");
      }
    } else {
      y_coords = read_coords(y_, ny);
    }

    // The ROI side decides which grid element type is used.
    if (mode_x == MODE_MESH) {
      MyMesh x_mesh;
      if (IOMesh<MyMesh>::vcgReadR(x_mesh, Rcpp::wrap(x_), x_it_) != 0) {
        Rcpp::stop("vcgDetectCollision: `x` has no faces; cannot use mesh mode.");
      }
      if (include_interior) {
        // Is_Inside needs per-face normals and the mesh bounding box.
        x_mesh.face.EnableNormal();
        vcg::tri::UpdateBounding<MyMesh>::Box(x_mesh);
        vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(x_mesh);
      } else {
        vcg::tri::UpdateBounding<MyMesh>::Box(x_mesh);
      }
      return dispatch_target<MyFace>(
        x_mesh.face.begin(), x_mesh.face.end(),
        static_cast<std::size_t>(x_mesh.fn), &x_mesh, radius,
        y_coords, ny, &y_mesh, mode_y, test_level, include_interior,
        x_mesh.fn);
    }

    int nx = 0;
    const std::vector<Scalar> x_coords = read_coords(x_, nx);

    if (mode_x == MODE_SEGMENTS) {
      // One grouping serves both the elements indexed and the chain count.
      const std::vector<Group> x_groups = find_groups(x_coords, nx);
      std::vector<SegElem> elems = build_seg_elems(x_coords, x_groups);
      return dispatch_target<SegElem>(
        elems.begin(), elems.end(), elems.size(), NULL, radius,
        y_coords, ny, &y_mesh, mode_y, test_level, include_interior,
        static_cast<int>(x_groups.size()));
    }

    std::vector<PointElem> elems = build_point_elems(x_coords, nx);
    return dispatch_target<PointElem>(
      elems.begin(), elems.end(), elems.size(), NULL, radius,
      y_coords, ny, &y_mesh, mode_y, test_level, include_interior, nx);

  } catch (std::exception &e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue;   // -Wall
}
