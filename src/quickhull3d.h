// Port of quickhull3d 3.1.2 to C++ (https://github.com/mauriciopoppe/quickhull3d,
// commit bef4fb4), for vcg_convex_hull(method = "quickhull3d").
//
// quickhull3d is Mauricio Poppe's port of John E. Lloyd's QuickHull3D
// (https://www.cs.ubc.ca/~lloyd/java/quickhull3d.html): Quickhull with faces
// as polygons in a half-edge structure, where neighboring faces that are
// coplanar or not convex within the tolerance are merged after every point
// added. The port keeps quickhull3d's arithmetic (the operation order of
// gl-vec3, get-plane-normal and point-line-distance) and its decisions, so it
// returns the same faces. It differs in that:
//   * the planar fallback (a 2D hull) is not ported; coplanar input throws;
//   * computeHorizon's recursion is an explicit stack that visits the faces
//     in the same order, so a large visible region cannot overflow the stack;
//   * the caller can raise the tolerance (not lower it);
//   * loops that a consistent run bounds (one point added per iteration, one
//     face deleted per merge) throw when they exceed the bound, instead of
//     running forever;
//   * two branches quickhull3d never reaches are left out: the minimum-area
//     normal correction (every face is created with minArea = 0) and the
//     NonConvex merge type (both merge passes use NonConvexWrtLargerFace,
//     whereas Lloyd's second pass uses NONCONVEX);
//   * JavaScript errors, including the TypeErrors that collinear or identical
//     points cause, become std::runtime_error.
//
// quickhull3d is distributed under the MIT license:
//
//   Copyright (c) Mauricio Poppe
//
//   Permission is hereby granted, free of charge, to any person obtaining
//   a copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//
//   The above copyright notice and this permission notice shall be
//   included in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
//   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef RAVETOOLS_QUICKHULL3D_H
#define RAVETOOLS_QUICKHULL3D_H

#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <stdexcept>
#include <vector>

// JavaScript rounds every multiplication and addition separately; clang by
// default fuses a * b + c into one rounding, which changes the last bits and,
// on inputs with exact ties, which faces come out
#if defined(__clang__)
#pragma clang fp contract(off)
#endif

namespace ravetools {
namespace quickhull3d {

// ---- gl-vec3, get-plane-normal and point-line-distance --------------------

inline double dot(const double *a, const double *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void add(double *out, const double *a, const double *b) {
  out[0] = a[0] + b[0];
  out[1] = a[1] + b[1];
  out[2] = a[2] + b[2];
}

inline void subtract(double *out, const double *a, const double *b) {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

inline void scale(double *out, const double *a, double s) {
  out[0] = a[0] * s;
  out[1] = a[1] * s;
  out[2] = a[2] * s;
}

inline void cross(double *out, const double *a, const double *b) {
  const double ax = a[0], ay = a[1], az = a[2];
  const double bx = b[0], by = b[1], bz = b[2];
  out[0] = ay * bz - az * by;
  out[1] = az * bx - ax * bz;
  out[2] = ax * by - ay * bx;
}

inline double squaredLength(const double *a) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

// leaves a zero vector unchanged
inline void normalize(double *out, const double *a) {
  double len = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
  if (len > 0) {
    len = 1 / std::sqrt(len);
    out[0] = a[0] * len;
    out[1] = a[1] * len;
    out[2] = a[2] * len;
  }
}

// get-plane-normal: normalize((p1 - p2) x (p2 - p3))
inline void planeNormal(double *out, const double *p1, const double *p2,
                        const double *p3) {
  double tmp[3];
  subtract(out, p1, p2);
  subtract(tmp, p2, p3);
  cross(out, out, tmp);
  normalize(out, out);
}

// point-line-distance: distance from p to the line through a and b
inline double pointLineDistance(const double *p, const double *a,
                                const double *b) {
  double ab[3], ap[3], cr[3];
  subtract(ab, b, a);
  subtract(ap, p, a);
  cross(cr, ap, ab);
  const double area = squaredLength(cr);
  const double s = squaredLength(ab);
  if (s == 0) {
    throw std::runtime_error("quickhull3d: all the points are identical");
  }
  return std::sqrt(area / s);
}

// ---- Vertex, VertexList, HalfEdge, Face ------------------------------------

struct Face;

struct Vertex {
  const double *point;  // x, y, z
  int index;            // position in the input
  Vertex *next;
  Vertex *prev;
  Face *face;           // the face that can see this point
  Vertex(const double *point, int index)
    : point(point), index(index), next(nullptr), prev(nullptr), face(nullptr) {}
};

// doubly linked list threaded through the vertices
class VertexList {
public:
  VertexList() : head(nullptr), tail(nullptr) {}

  void clear() { head = tail = nullptr; }

  void insertBefore(Vertex *target, Vertex *node) {
    node->prev = target->prev;
    node->next = target;
    if (!node->prev) {
      head = node;
    } else {
      node->prev->next = node;
    }
    target->prev = node;
  }

  void add(Vertex *node) {
    if (!head) {
      head = node;
    } else {
      tail->next = node;
    }
    node->prev = tail;
    node->next = nullptr;
    tail = node;
  }

  // appends the chain that starts at node
  void addAll(Vertex *node) {
    if (!head) {
      head = node;
    } else {
      tail->next = node;
    }
    node->prev = tail;
    while (node->next) {
      node = node->next;
    }
    tail = node;
  }

  void remove(Vertex *node) {
    if (!node->prev) {
      head = node->next;
    } else {
      node->prev->next = node->next;
    }
    if (!node->next) {
      tail = node->prev;
    } else {
      node->next->prev = node->prev;
    }
  }

  // removes the chain from a to b, a before b
  void removeChain(Vertex *a, Vertex *b) {
    if (!a->prev) {
      head = b->next;
    } else {
      a->prev->next = b->next;
    }
    if (!b->next) {
      tail = a->prev;
    } else {
      b->next->prev = a->prev;
    }
  }

  Vertex *first() const { return head; }
  bool isEmpty() const { return !head; }

private:
  Vertex *head;
  Vertex *tail;
};

struct HalfEdge {
  Vertex *vertex;  // the head
  Face *face;
  HalfEdge *next;
  HalfEdge *prev;
  HalfEdge *opposite;

  HalfEdge(Vertex *vertex, Face *face)
    : vertex(vertex), face(face), next(nullptr), prev(nullptr),
      opposite(nullptr) {}

  Vertex *head() const { return vertex; }
  Vertex *tail() const { return prev ? prev->vertex : nullptr; }

  void setOpposite(HalfEdge *edge) {
    opposite = edge;
    edge->opposite = this;
  }
};

enum Mark { Visible = 0, NonConvex = 1, Deleted = 2 };

struct Face {
  double normal[3];
  double centroid[3];
  double offset;     // signed distance from the origin to the plane
  Vertex *outside;   // first of the claimed vertices this face can see
  Mark mark;
  HalfEdge *edge;
  int nVertices;
  double area;

  Face()
    : normal{0, 0, 0}, centroid{0, 0, 0}, offset(0), outside(nullptr),
      mark(Visible), edge(nullptr), nVertices(0), area(0) {}

  HalfEdge *getEdge(int i) const {
    HalfEdge *it = edge;
    while (i > 0) {
      it = it->next;
      i -= 1;
    }
    while (i < 0) {
      it = it->prev;
      i += 1;
    }
    return it;
  }

  // Newell-style sum over the fan from the first vertex; a zero area gives a
  // NaN normal, as in quickhull3d
  void computeNormal() {
    HalfEdge *e0 = edge;
    HalfEdge *e1 = e0->next;
    HalfEdge *e2 = e1->next;
    double v1[3], v2[3], t[3];
    subtract(v2, e1->head()->point, e0->head()->point);
    nVertices = 2;
    normal[0] = normal[1] = normal[2] = 0;
    while (e2 != e0) {
      v1[0] = v2[0];
      v1[1] = v2[1];
      v1[2] = v2[2];
      subtract(v2, e2->head()->point, e0->head()->point);
      cross(t, v1, v2);
      add(normal, normal, t);
      e2 = e2->next;
      nVertices += 1;
    }
    area = std::sqrt(squaredLength(normal));
    scale(normal, normal, 1 / area);
  }

  void computeCentroid() {
    centroid[0] = centroid[1] = centroid[2] = 0;
    HalfEdge *e = edge;
    do {
      add(centroid, centroid, e->head()->point);
      e = e->next;
    } while (e != edge);
    scale(centroid, centroid, 1.0 / nVertices);
  }

  void computeNormalAndCentroid() {
    computeNormal();
    computeCentroid();
    offset = dot(normal, centroid);
  }

  double distanceToPlane(const double *point) const {
    return dot(normal, point) - offset;
  }

  // Connects prev to next (prev's head is next's tail). When both border the
  // same face, the edge between them is redundant and is removed, and a
  // triangular neighbor that loses an edge that way is returned to discard.
  Face *connectHalfEdges(HalfEdge *prev, HalfEdge *next) {
    Face *discardedFace = nullptr;
    if (prev->opposite->face == next->opposite->face) {
      Face *oppositeFace = next->opposite->face;
      HalfEdge *oppositeEdge;
      if (prev == edge) {
        edge = next;
      }
      if (oppositeFace->nVertices == 3) {
        oppositeEdge = next->opposite->prev->opposite;
        oppositeFace->mark = Deleted;
        discardedFace = oppositeFace;
      } else {
        oppositeEdge = next->opposite->next;
        if (oppositeFace->edge == oppositeEdge->prev) {
          oppositeFace->edge = oppositeEdge;
        }
        oppositeEdge->prev = oppositeEdge->prev->prev;
        oppositeEdge->prev->next = oppositeEdge;
      }
      next->prev = prev->prev;
      next->prev->next = next;
      next->setOpposite(oppositeEdge);
      oppositeFace->computeNormalAndCentroid();
    } else {
      prev->next = next;
      next->prev = prev;
    }
    return discardedFace;
  }

  // Absorbs the face across adjacentEdge; appends the faces destroyed on the
  // way to discardedFaces
  void mergeAdjacentFaces(HalfEdge *adjacentEdge,
                          std::vector<Face *> &discardedFaces) {
    HalfEdge *oppositeEdge = adjacentEdge->opposite;
    Face *oppositeFace = oppositeEdge->face;

    discardedFaces.push_back(oppositeFace);
    oppositeFace->mark = Deleted;

    // the chain of edges shared with oppositeFace
    HalfEdge *adjacentEdgePrev = adjacentEdge->prev;
    HalfEdge *adjacentEdgeNext = adjacentEdge->next;
    HalfEdge *oppositeEdgePrev = oppositeEdge->prev;
    HalfEdge *oppositeEdgeNext = oppositeEdge->next;
    while (adjacentEdgePrev->opposite->face == oppositeFace) {
      adjacentEdgePrev = adjacentEdgePrev->prev;
      oppositeEdgeNext = oppositeEdgeNext->next;
    }
    while (adjacentEdgeNext->opposite->face == oppositeFace) {
      adjacentEdgeNext = adjacentEdgeNext->next;
      oppositeEdgePrev = oppositeEdgePrev->prev;
    }

    // the edges of oppositeFace that are not shared now belong to this face
    for (HalfEdge *e = oppositeEdgeNext; e != oppositeEdgePrev->next; e = e->next) {
      e->face = this;
    }

    // a next edge, because connectHalfEdges can destroy prev edges
    edge = adjacentEdgeNext;

    Face *discardedFace = connectHalfEdges(oppositeEdgePrev, adjacentEdgeNext);
    if (discardedFace) {
      discardedFaces.push_back(discardedFace);
    }
    discardedFace = connectHalfEdges(adjacentEdgePrev, oppositeEdgeNext);
    if (discardedFace) {
      discardedFaces.push_back(discardedFace);
    }

    computeNormalAndCentroid();
  }

  void collectIndices(std::vector<int> &indices) const {
    indices.clear();
    HalfEdge *e = edge;
    do {
      indices.push_back(e->head()->index);
      e = e->next;
    } while (e != edge);
  }
};

// ---- QuickHull --------------------------------------------------------------

class QuickHull {
public:
  // points: n consecutive (x, y, z) triplets, which must outlive the object;
  // a negative tolerance uses quickhull3d's own
  QuickHull(const double *points, std::size_t n, double tolerance = -1)
    : tolerance(-1), toleranceOverride(tolerance) {
    if (n < 4) {
      throw std::runtime_error("quickhull3d: cannot build a simplex out of <4 points");
    }
    vertices.reserve(n);
    for (std::size_t i = 0; i < n; i += 1) {
      vertices.push_back(Vertex(points + 3 * i, (int)i));
    }
  }

  void build() {
    Vertex *v[4];
    computeTetrahedronExtremes(v);
    if (allPointsBelongToPlane(v[0], v[1], v[2])) {
      throw std::runtime_error("quickhull3d: the points are coplanar");
    }
    createInitialSimplex(v[0], v[1], v[2], v[3]);
    Vertex *eyeVertex;
    std::size_t iterations = 0;
    while ((eyeVertex = nextVertexToAdd())) {
      // every point is added at most once
      if (++iterations > vertices.size()) {
        throw std::runtime_error("quickhull3d: the hull did not converge");
      }
      addVertexToHull(eyeVertex);
    }
    reindexFaceAndVertices();
  }

  // the faces, each split into a fan of triangles from its first vertex
  std::vector<std::array<int, 3> > collectFaces() const {
    std::vector<std::array<int, 3> > triangles;
    std::vector<int> indices;
    for (std::size_t i = 0; i < faces.size(); i += 1) {
      if (faces[i]->mark != Visible) {
        throw std::runtime_error("quickhull3d: attempt to include a destroyed face in the hull");
      }
      faces[i]->collectIndices(indices);
      for (std::size_t j = 0; j + 2 < indices.size(); j += 1) {
        std::array<int, 3> triangle = {{ indices[0], indices[j + 1], indices[j + 2] }};
        triangles.push_back(triangle);
      }
    }
    return triangles;
  }

  double tolerance;

private:
  double toleranceOverride;
  std::vector<Vertex> vertices;
  std::deque<Face> facePool;      // stable addresses; nothing is freed early
  std::deque<HalfEdge> edgePool;
  std::vector<Face *> faces;
  std::vector<Face *> newFaces;
  VertexList claimed;
  VertexList unclaimed;

  Face *createTriangle(Vertex *v0, Vertex *v1, Vertex *v2) {
    facePool.emplace_back();
    Face *face = &facePool.back();
    edgePool.emplace_back(v0, face);
    HalfEdge *e0 = &edgePool.back();
    edgePool.emplace_back(v1, face);
    HalfEdge *e1 = &edgePool.back();
    edgePool.emplace_back(v2, face);
    HalfEdge *e2 = &edgePool.back();
    e0->next = e2->prev = e1;
    e1->next = e0->prev = e2;
    e2->next = e1->prev = e0;
    face->edge = e0;
    face->computeNormalAndCentroid();
    return face;
  }

  void addVertexToFace(Vertex *vertex, Face *face) {
    vertex->face = face;
    if (!face->outside) {
      claimed.add(vertex);
    } else {
      claimed.insertBefore(face->outside, vertex);
    }
    face->outside = vertex;
  }

  void removeVertexFromFace(Vertex *vertex, Face *face) {
    if (vertex == face->outside) {
      if (vertex->next && vertex->next->face == face) {
        face->outside = vertex->next;
      } else {
        face->outside = nullptr;
      }
    }
    claimed.remove(vertex);
  }

  // detaches the claimed vertices that face can see; returns the first
  Vertex *removeAllVerticesFromFace(Face *face) {
    if (face->outside) {
      Vertex *end = face->outside;
      while (end->next && end->next->face == face) {
        end = end->next;
      }
      claimed.removeChain(face->outside, end);
      end->next = nullptr;
      return face->outside;
    }
    return nullptr;
  }

  // moves the vertices face can see to absorbingFace when they can see it,
  // otherwise to the unclaimed list
  void deleteFaceVertices(Face *face, Face *absorbingFace) {
    Vertex *faceVertices = removeAllVerticesFromFace(face);
    if (faceVertices) {
      if (!absorbingFace) {
        unclaimed.addAll(faceVertices);
      } else {
        Vertex *nextVertex;
        for (Vertex *vertex = faceVertices; vertex; vertex = nextVertex) {
          nextVertex = vertex->next;
          const double distance = absorbingFace->distanceToPlane(vertex->point);
          if (distance > tolerance) {
            addVertexToFace(vertex, absorbingFace);
          } else {
            unclaimed.add(vertex);
          }
        }
      }
    }
  }

  void resolveUnclaimedPoints(const std::vector<Face *> &candidates) {
    Vertex *vertexNext = unclaimed.first();
    for (Vertex *vertex = vertexNext; vertex; vertex = vertexNext) {
      vertexNext = vertex->next;
      double maxDistance = tolerance;
      Face *maxFace = nullptr;
      for (std::size_t i = 0; i < candidates.size(); i += 1) {
        Face *face = candidates[i];
        if (face->mark == Visible) {
          const double dist = face->distanceToPlane(vertex->point);
          if (dist > maxDistance) {
            maxDistance = dist;
            maxFace = face;
          }
          if (maxDistance > 1000 * tolerance) {
            break;
          }
        }
      }
      if (maxFace) {
        addVertexToFace(vertex, maxFace);
      }
    }
  }

  bool allPointsBelongToPlane(Vertex *v0, Vertex *v1, Vertex *v2) const {
    double normal[3];
    planeNormal(normal, v0->point, v1->point, v2->point);
    const double distToPlane = dot(normal, v0->point);
    for (std::size_t i = 0; i < vertices.size(); i += 1) {
      const double dist = dot(vertices[i].point, normal);
      if (std::fabs(dist - distToPlane) > tolerance) {
        return false;
      }
    }
    return true;
  }

  // the extremes of the initial tetrahedron; also sets the tolerance
  void computeTetrahedronExtremes(Vertex **v) {
    double min[3], max[3];
    Vertex *minVertices[3];
    Vertex *maxVertices[3];
    for (int i = 0; i < 3; i += 1) {
      minVertices[i] = maxVertices[i] = &vertices[0];
      min[i] = max[i] = vertices[0].point[i];
    }
    for (std::size_t i = 1; i < vertices.size(); i += 1) {
      Vertex *vertex = &vertices[i];
      const double *point = vertex->point;
      for (int j = 0; j < 3; j += 1) {
        if (point[j] < min[j]) {
          min[j] = point[j];
          minVertices[j] = vertex;
        }
      }
      for (int j = 0; j < 3; j += 1) {
        if (point[j] > max[j]) {
          max[j] = point[j];
          maxVertices[j] = vertex;
        }
      }
    }

    tolerance =
      3 * std::numeric_limits<double>::epsilon() *
      (std::fmax(std::fabs(min[0]), std::fabs(max[0])) +
       std::fmax(std::fabs(min[1]), std::fabs(max[1])) +
       std::fmax(std::fabs(min[2]), std::fabs(max[2])));
    // A caller can only raise this tolerance: below the rounding error of the
    // coordinates the merge decisions contradict each other and the merge
    // loop never ends (in quickhull3d itself too, with its tolerance lowered)
    if (toleranceOverride > tolerance) {
      tolerance = toleranceOverride;
    }

    // the two vertices with the greatest separation along one axis
    double maxDistance = 0;
    int indexMax = 0;
    for (int i = 0; i < 3; i += 1) {
      const double distance = maxVertices[i]->point[i] - minVertices[i]->point[i];
      if (distance > maxDistance) {
        maxDistance = distance;
        indexMax = i;
      }
    }
    Vertex *v0 = minVertices[indexMax];
    Vertex *v1 = maxVertices[indexMax];
    Vertex *v2 = nullptr;
    Vertex *v3 = nullptr;

    // the vertex farthest from the line through v0 and v1
    maxDistance = 0;
    for (std::size_t i = 0; i < vertices.size(); i += 1) {
      Vertex *vertex = &vertices[i];
      if (vertex != v0 && vertex != v1) {
        const double distance = pointLineDistance(vertex->point, v0->point, v1->point);
        if (distance > maxDistance) {
          maxDistance = distance;
          v2 = vertex;
        }
      }
    }
    if (!v2) {
      throw std::runtime_error("quickhull3d: the points are collinear");
    }

    // the vertex farthest from the plane through v0, v1 and v2
    double normal[3];
    planeNormal(normal, v0->point, v1->point, v2->point);
    const double distPO = dot(v0->point, normal);
    maxDistance = -1;
    for (std::size_t i = 0; i < vertices.size(); i += 1) {
      Vertex *vertex = &vertices[i];
      if (vertex != v0 && vertex != v1 && vertex != v2) {
        const double distance = std::fabs(dot(normal, vertex->point) - distPO);
        if (distance > maxDistance) {
          maxDistance = distance;
          v3 = vertex;
        }
      }
    }

    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
  }

  // the initial tetrahedron, with every other point assigned to the face it
  // is farthest above
  void createInitialSimplex(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3) {
    double normal[3];
    planeNormal(normal, v0->point, v1->point, v2->point);
    const double distPO = dot(v0->point, normal);

    Face *simplex[4];
    if (dot(v3->point, normal) - distPO < 0) {
      // the plane normal points outside the tetrahedron
      simplex[0] = createTriangle(v0, v1, v2);
      simplex[1] = createTriangle(v3, v1, v0);
      simplex[2] = createTriangle(v3, v2, v1);
      simplex[3] = createTriangle(v3, v0, v2);
      for (int i = 0; i < 3; i += 1) {
        const int j = (i + 1) % 3;
        simplex[i + 1]->getEdge(2)->setOpposite(simplex[0]->getEdge(j));
        simplex[i + 1]->getEdge(1)->setOpposite(simplex[j + 1]->getEdge(0));
      }
    } else {
      // the plane normal points inside the tetrahedron
      simplex[0] = createTriangle(v0, v2, v1);
      simplex[1] = createTriangle(v3, v0, v1);
      simplex[2] = createTriangle(v3, v1, v2);
      simplex[3] = createTriangle(v3, v2, v0);
      for (int i = 0; i < 3; i += 1) {
        const int j = (i + 1) % 3;
        simplex[i + 1]->getEdge(2)->setOpposite(simplex[0]->getEdge((3 - i) % 3));
        simplex[i + 1]->getEdge(0)->setOpposite(simplex[j + 1]->getEdge(1));
      }
    }

    for (int i = 0; i < 4; i += 1) {
      faces.push_back(simplex[i]);
    }

    for (std::size_t i = 0; i < vertices.size(); i += 1) {
      Vertex *vertex = &vertices[i];
      if (vertex != v0 && vertex != v1 && vertex != v2 && vertex != v3) {
        double maxDistance = tolerance;
        Face *maxFace = nullptr;
        for (int j = 0; j < 4; j += 1) {
          const double distance = simplex[j]->distanceToPlane(vertex->point);
          if (distance > maxDistance) {
            maxDistance = distance;
            maxFace = simplex[j];
          }
        }
        if (maxFace) {
          addVertexToFace(vertex, maxFace);
        }
      }
    }
  }

  void reindexFaceAndVertices() {
    std::vector<Face *> activeFaces;
    for (std::size_t i = 0; i < faces.size(); i += 1) {
      if (faces[i]->mark == Visible) {
        activeFaces.push_back(faces[i]);
      }
    }
    faces.swap(activeFaces);
  }

  // the farthest vertex seen by the face of the first claimed vertex
  Vertex *nextVertexToAdd() {
    if (claimed.isEmpty()) {
      return nullptr;
    }
    Vertex *eyeVertex = nullptr;
    double maxDistance = 0;
    Face *eyeFace = claimed.first()->face;
    for (Vertex *vertex = eyeFace->outside; vertex && vertex->face == eyeFace;
         vertex = vertex->next) {
      const double distance = eyeFace->distanceToPlane(vertex->point);
      if (distance > maxDistance) {
        maxDistance = distance;
        eyeVertex = vertex;
      }
    }
    return eyeVertex;
  }

  // Deletes the faces that can see eyePoint, starting from face, and collects
  // in horizon, counter-clockwise, the edges between a deleted face and one
  // that stays. quickhull3d recurses into each face across an edge; the
  // stack holds, for every face entered, the edge to visit next and the edge
  // it was entered through (where its visit ends).
  void computeHorizon(const double *eyePoint, Face *face,
                      std::vector<HalfEdge *> &horizon) {
    struct Frame {
      HalfEdge *edge;
      HalfEdge *crossEdge;
    };
    std::vector<Frame> stack;

    deleteFaceVertices(face, nullptr);
    face->mark = Deleted;
    HalfEdge *first = face->getEdge(0);
    stack.push_back(Frame{ first, first });

    while (!stack.empty()) {
      HalfEdge *edge = stack.back().edge;
      HalfEdge *next = edge->next;
      const bool last = next == stack.back().crossEdge;
      stack.back().edge = next;
      if (last) {
        stack.pop_back();
      }

      HalfEdge *oppositeEdge = edge->opposite;
      Face *oppositeFace = oppositeEdge->face;
      if (oppositeFace->mark == Visible) {
        if (oppositeFace->distanceToPlane(eyePoint) > tolerance) {
          deleteFaceVertices(oppositeFace, nullptr);
          oppositeFace->mark = Deleted;
          stack.push_back(Frame{ oppositeEdge->next, oppositeEdge });
        } else {
          horizon.push_back(edge);
        }
      }
    }
  }

  // the triangle (eyeVertex, horizonEdge's tail, its head), joined to the
  // face across the horizon; returns its edge whose head is eyeVertex
  HalfEdge *addAdjoiningFace(Vertex *eyeVertex, HalfEdge *horizonEdge) {
    Face *face = createTriangle(eyeVertex, horizonEdge->tail(), horizonEdge->head());
    faces.push_back(face);
    face->getEdge(-1)->setOpposite(horizonEdge->opposite);
    return face->getEdge(0);
  }

  void addNewFaces(Vertex *eyeVertex, const std::vector<HalfEdge *> &horizon) {
    if (horizon.empty()) {
      throw std::runtime_error("quickhull3d: empty horizon");
    }
    newFaces.clear();
    HalfEdge *firstSideEdge = nullptr;
    HalfEdge *previousSideEdge = nullptr;
    for (std::size_t i = 0; i < horizon.size(); i += 1) {
      HalfEdge *sideEdge = addAdjoiningFace(eyeVertex, horizon[i]);
      if (!firstSideEdge) {
        firstSideEdge = sideEdge;
      } else {
        sideEdge->next->setOpposite(previousSideEdge);
      }
      newFaces.push_back(sideEdge->face);
      previousSideEdge = sideEdge;
    }
    firstSideEdge->next->setOpposite(previousSideEdge);
  }

  // distance from the centroid of the face across edge to edge's face
  // (positive when the two are concave)
  double oppositeFaceDistance(HalfEdge *edge) const {
    return edge->face->distanceToPlane(edge->opposite->face->centroid);
  }

  // Merges face with the first neighbor that is concave or coplanar with it
  // within the tolerance, judged from the larger face; marks face NonConvex
  // when only the smaller face's view says they are not convex
  bool doAdjacentMerge(Face *face) {
    HalfEdge *edge = face->edge;
    bool convex = true;
    int it = 0;
    do {
      if (it >= face->nVertices) {
        throw std::runtime_error("quickhull3d: merge recursion limit exceeded");
      }
      Face *oppositeFace = edge->opposite->face;
      bool merge = false;
      if (face->area > oppositeFace->area) {
        if (oppositeFaceDistance(edge) > -tolerance) {
          merge = true;
        } else if (oppositeFaceDistance(edge->opposite) > -tolerance) {
          convex = false;
        }
      } else {
        if (oppositeFaceDistance(edge->opposite) > -tolerance) {
          merge = true;
        } else if (oppositeFaceDistance(edge) > -tolerance) {
          convex = false;
        }
      }

      if (merge) {
        std::vector<Face *> discardedFaces;
        face->mergeAdjacentFaces(edge, discardedFaces);
        for (std::size_t i = 0; i < discardedFaces.size(); i += 1) {
          deleteFaceVertices(discardedFaces[i], face);
        }
        return true;
      }

      edge = edge->next;
      it += 1;
    } while (edge != face->edge);
    if (!convex) {
      face->mark = NonConvex;
    }
    return false;
  }

  void addVertexToHull(Vertex *eyeVertex) {
    std::vector<HalfEdge *> horizon;
    unclaimed.clear();

    // so that eyeVertex does not land in the unclaimed list
    removeVertexFromFace(eyeVertex, eyeVertex->face);
    computeHorizon(eyeVertex->point, eyeVertex->face, horizon);
    addNewFaces(eyeVertex, horizon);

    // every merge deletes a face, so a consistent run merges fewer times than
    // there are faces
    std::size_t merges = 0;
    const std::size_t maxMerges = faces.size();

    // first pass: merge with respect to the larger face
    for (std::size_t i = 0; i < newFaces.size(); i += 1) {
      Face *face = newFaces[i];
      if (face->mark == Visible) {
        while (doAdjacentMerge(face)) {
          if (++merges > maxMerges) {
            throw std::runtime_error("quickhull3d: face merging did not terminate");
          }
        }
      }
    }
    // second pass: faces the first pass marked as not convex
    for (std::size_t i = 0; i < newFaces.size(); i += 1) {
      Face *face = newFaces[i];
      if (face->mark == NonConvex) {
        face->mark = Visible;
        while (doAdjacentMerge(face)) {
          if (++merges > maxMerges) {
            throw std::runtime_error("quickhull3d: face merging did not terminate");
          }
        }
      }
    }

    resolveUnclaimedPoints(newFaces);
  }
};

} // namespace quickhull3d
} // namespace ravetools

#endif // RAVETOOLS_QUICKHULL3D_H
