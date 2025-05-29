#ifndef RAVETOOLS_VCG_COMMON_H
#define RAVETOOLS_VCG_COMMON_H

#ifndef NDEBUG
#define NDEBUG
#endif

#include <string.h>
#include <vector>
#include <stdio.h>
#include <cstddef>
#include <stdexcept>
#include <Rcpp.h>

#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/update/topology.h>

//#include <vcg/complex/algorithms/update/edges.h>

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/space/point3.h>

#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/container/simple_temporary_data.h>

#include <vcg/complex/algorithms/intersection.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/resampler.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>
#include<vcg/complex/algorithms/create/platonic.h>
//include headers for search grids
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/implicit_smooth.h>
#include <vcg/complex/algorithms/geodesic.h>
//#include <vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>
#include <vcg/container/simple_temporary_data.h>

// VCG File Format Importer/Exporter
// #include <wrap/io_trimesh/import.h>
// #include <wrap/io_trimesh/import_off.h>
// #include <wrap/io_trimesh/export.h>
// #include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/refine_loop.h>

//#include <vcg/complex/algorithms/update/curvature.h>


namespace ravetools
{

class MyFace;
class MyEdge;
class MyVertex;



struct MyUsedTypes: public vcg::UsedTypes<
  vcg::Use<MyVertex>::AsVertexType,
  vcg::Use<MyEdge>::AsEdgeType,
  vcg::Use<MyFace>::AsFaceType
> {};

class MyEdge : public vcg::Edge<MyUsedTypes> {};

class MyVertex : public vcg::Vertex<
  MyUsedTypes,
  vcg::vertex::InfoOcf,
  vcg::vertex::Coord3f,
  vcg::vertex::BitFlags,
  vcg::vertex::Normal3f,
  vcg::vertex::Mark,
  vcg::vertex::Color4bOcf,
  vcg::vertex::QualityfOcf,
  vcg::vertex::VFAdjOcf,
  // vcg::vertex::CurvatureOcf,
  vcg::vertex::CurvatureDirfOcf,
  vcg::vertex::TexCoordfOcf
 /* vcg::vertex::VFAdj,
  vcg::vertex::Curvaturef,
  vcg::vertex::CurvatureDirf */
> {};

class MyFace: public vcg::Face <
  MyUsedTypes,
  vcg::face::InfoOcf,
  vcg::face::VertexRef,
  vcg::face::BitFlags,
  vcg::face::Mark,
  vcg::face::FFAdjOcf,
  vcg::face::VFAdjOcf,
  vcg::face::WedgeTexCoordfOcf,
  vcg::face::Color4bOcf,
  vcg::face::QualityfOcf,
  vcg::face::Normal3fOcf
> {};


// ocf class
class MyMesh : public vcg::tri::TriMesh<
  vcg::vertex::vector_ocf<MyVertex>,
  vcg::face::vector_ocf<MyFace >
> {};

// default class
//class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};
typedef  MyMesh::ScalarType ScalarType;
typedef  MyMesh::VertexIterator VertexIterator;
typedef  MyMesh::VertexPointer VertexPointer;
typedef  MyMesh::FaceIterator FaceIterator;
typedef  MyMesh::FacePointer FacePointer;
typedef  MyMesh::EdgePointer   EdgePointer;
typedef  MyMesh::EdgeIterator   EdgeIterator;
typedef  MyMesh::CoordType CoordType;
typedef  MyMesh::ScalarType ScalarType;
typedef  MyMesh::VertContainer VertContainer;
typedef  MyMesh::FaceContainer FaceContainer;
/*typedef MyMesh::ConstVertexIterator ConstVertexIterator;
typedef MyMesh::ConstFaceIterator   ConstFaceIterator;*/

// for point clouds
class MyPointCloudVertex;
struct MyPointCloudUsedTypes : public vcg::UsedTypes<
  vcg::Use<MyPointCloudVertex>::AsVertexType
>{};

class MyPointCloudEdge : public vcg::Edge<MyPointCloudUsedTypes> {};

class MyPointCloudVertex : public vcg::Vertex<
  MyPointCloudUsedTypes,
  // vertex::InfoOcf,
  vcg::vertex::Coord3f,
  vcg::vertex::BitFlags,
  vcg::vertex::Normal3f,
  vcg::vertex::Mark,
  vcg::vertex::Color4b,
  vcg::vertex::Qualityf
>{};

class MyPointCloud : public vcg::tri::TriMesh< std::vector<MyPointCloudVertex> > {};



template <class VOX_TYPE>
class MySimpleVolume : public vcg::BasicGrid<float>
{
public:
  typedef VOX_TYPE VoxelType;

  const vcg::Point3i &ISize() { return siz; }

  float Val( const int &x, const int &y, const int &z ) const {
    return cV(x,y,z).V();
  }

  float &Val(const int &x,const int &y,const int &z) {
    return V(x,y,z).V();
    //else return numeric_limits<float>::quiet_NaN( );
  }

  VOX_TYPE &V(const int &x,const int &y,const int &z) {
    return Vol[x+y*siz[0]+z*siz[0]*siz[1]];
  }

  VOX_TYPE &V(const vcg::Point3i &pi) {
    return Vol[ pi[0] + pi[1]*siz[0] + pi[2]*siz[0]*siz[1] ];
  }

  const VOX_TYPE &cV(const int &x,const int &y,const int &z) const {
    return Vol[x+y*siz[0]+z*siz[0]*siz[1]];
  }
  bool ValidCell(const vcg::Point3i & /*p0*/, const vcg::Point3i & /*p1*/) const { return true;}

  template < class VertexPointerType >
  void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,XAxis>(p1,p2,v,thr); }

  template < class VertexPointerType >
  void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,YAxis>(p1,p2,v,thr); }

  template < class VertexPointerType >
  void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,ZAxis>(p1,p2,v,thr); }

  /// The following members/methods are just for this particular case.
  /// The above one are the one required by the marching cube interface.

  std::vector<VoxelType> Vol;

  typedef enum { XAxis=0,YAxis=1,ZAxis=2} VolumeAxis;

  template < class VertexPointerType,  VolumeAxis AxisVal >
  void GetIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  {
    float f1 = V(p1).V()-thr;
    float f2 = V(p2).V()-thr;
    float u = (float) f1/(f1-f2);
    if(AxisVal==XAxis) v->P().X() = (float) p1.X()*(1-u) + u*p2.X();
    else v->P().X() = (float) p1.X();
    if(AxisVal==YAxis) v->P().Y() = (float) p1.Y()*(1-u) + u*p2.Y();
    else v->P().Y() = (float) p1.Y();
    if(AxisVal==ZAxis) v->P().Z() = (float) p1.Z()*(1-u) + u*p2.Z();
    else v->P().Z() = (float) p1.Z();
    //this->IPfToPf(v->P(),v->P());
    if(VoxelType::HasNormal()) v->N() = V(p1).N()*(1-u) + V(p2).N()*u;
  }



  void Init(vcg::Point3i _sz)
  {
    siz=_sz;
    Vol.resize(siz[0]*siz[1]*siz[2]);
  }



};

class MySimpleVoxel
{
private:
  float _v;
public:
  float &V() {return _v;}
  float V() const {return _v;}
  static bool HasNormal() {return false;}
  vcg::Point3f N() const {return vcg::Point3f(0,0,0);}
  vcg::Point3f &N()  { static vcg::Point3f _p(0,0,0); return _p;}
};


typedef MySimpleVolume<MySimpleVoxel> MyVolume;
typedef vcg::tri::TrivialWalker<MyMesh, MyVolume>	MyWalker;
typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;



template <class IOMeshType>
class IOMesh
{
public:
  typedef IOMeshType MeshType;
  typedef typename MeshType::CoordType      CoordType;
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef typename MeshType::FaceContainer  FaceContainer;
  typedef typename MeshType::VertContainer  VertContainer;

  // Fill an empty mesh with vertices and faces from R
  static int vcgReadR(
      MeshType &m, SEXP vb_, SEXP it_ = Rcpp::wrap(0), SEXP normals_ = Rcpp::wrap(0),
      bool zerobegin = true, bool readnormals = true, bool readfaces = true
  ) {
    try {
      //insert vertices
      if ( Rf_isMatrix(vb_) && VertexType::HasCoord() ) {
        Rcpp::NumericMatrix vb(vb_);
        int d =  vb.ncol();
        vcg::tri::Allocator<MeshType>::AddVertices(m,d);
        std::vector<VertexPointer> ivp;
        ivp.resize(d);
        vcg::SimpleTempData<typename MeshType::VertContainer, unsigned int> indices(m.vert);
        //read vertices
        for ( int i = 0; i < d; i++ ) {
          VertexIterator vi = m.vert.begin() + i;
          ivp[i] = &*vi;
          (*vi).P() = CoordType( vb(0,i), vb(1,i), vb(2,i) );
        }
        //insert vertex normals
        if ( Rf_isMatrix(normals_) && vcg::tri::HasPerVertexNormal(m) && readnormals ) {
          Rcpp::NumericMatrix normals(normals_);
          if ( normals.ncol() != d ) {
            Rprintf("number of normals is not equal to number of vertices");
          } else {
            vcg::SimpleTempData<typename MeshType::VertContainer, unsigned int> indices(m.vert);

            // #pragma omp parallel for schedule(static)
            for ( int i = 0; i < d; i++ ) {
              VertexIterator vi = m.vert.begin()+i;
              ivp[i] = &*vi;
              (*vi).N() = CoordType( normals(0,i), normals(1,i), normals(2,i) );
            }
          }
        }
        //process faces but check attributes and input first
        if ( Rf_isMatrix(it_) && FaceType::HasVertexRef() && readfaces ) {
          Rcpp::IntegerMatrix it(it_);
          unsigned int faced = it.ncol();
          vcg::tri::Allocator<MeshType>::AddFaces(m, faced);
          vcg::SimpleTempData<typename MeshType::FaceContainer, unsigned int> indicesf(m.face);

          for (unsigned int i=0; i < faced ; i++) {
            int subtract = 0;
            if (!zerobegin) {
              subtract = 1;
            }
            FaceIterator fi = m.face.begin()+i;
            indicesf[fi] = i;
            for (int j = 0; j < 3; j++)
              (*fi).V(j)=ivp[it(j,i)-subtract];
          }
          return 0;
        } else {
          return 1;
        }
      } else {
        return -1;
      }
    } catch (std::exception& e) {
      Rcpp::stop( e.what());
      return 1;
    } catch (...) {
      Rcpp::stop("unknown exception");
      return 1;
    }
  };

  static Rcpp::List vcgToR(MeshType &mesh, bool includeNormals = false) {
    try {

      Rcpp::List out;

      // Map each valid vertex iterator -> compact index
      vcg::SimpleTempData<typename MeshType::VertContainer,unsigned int> vertIndexMap(mesh.vert);

      // Preallocate full-size buffers
      Rcpp::NumericMatrix vertexMat(4, mesh.vn), normalMat(4, mesh.vn);
      std::fill(vertexMat.begin(), vertexMat.end(), 1.0);
      std::fill(normalMat.begin(), normalMat.end(), 1.0);

      // Preallocate triangle index buffer
      Rcpp::IntegerMatrix faceIdxMat(3, mesh.fn);

      // 1) Gather vertices
      VertexIterator vertIt;
      VertexPointer vertPtr;
      size_t validVertCount = 0;
      for (size_t i = 0; i < mesh.vn; ++i) {
        vertIt  = mesh.vert.begin() + i;
        vertPtr = &(*vertIt);
        if (vertPtr && !vertPtr->IsD()) {
          // assign a compact index
          vertIndexMap[vertIt] = validVertCount;

          // copy position...
          for (int coord = 0; coord < 3; ++coord) {
            vertexMat(coord, validVertCount) = vertPtr->P()[coord];
            if (includeNormals) {
              // ...and normal if requested
              normalMat(coord, validVertCount) = vertPtr->N()[coord];
            }
          }
          ++validVertCount;
        }
      }

      // 2) Gather faces
      FaceIterator faceIt;
      FacePointer facePtr;
      size_t validFaceCount = 0;
      for (size_t i = 0; i < mesh.fn; ++i) {
        faceIt  = mesh.face.begin() + i;
        facePtr = &(*faceIt);
        if (facePtr && !facePtr->IsD()
              && facePtr->V(0) && facePtr->V(1) && facePtr->V(2)) {
          // copy vertex indices (+1 for R’s 1-based)
          for (int corner = 0; corner < 3; ++corner) {
            faceIdxMat(corner, validFaceCount) = vertIndexMap[facePtr->cV(corner)] + 1;
          }
          ++validFaceCount;
        }
      }

      // 3) Slice off only the filled-in columns
      Rcpp::NumericMatrix vertexMatSub =
        vertexMat(      Rcpp::_, Rcpp::Range(0, validVertCount - 1));
      Rcpp::IntegerMatrix faceIdxSub =
        faceIdxMat(     Rcpp::_, Rcpp::Range(0, validFaceCount - 1));

      out["vb"] = vertexMatSub;
      out["it"] = faceIdxSub;

      if (includeNormals) {
        Rcpp::NumericMatrix normalMatSub =
          normalMat(Rcpp::_, Rcpp::Range(0, validVertCount - 1));
        out["normals"] = normalMatSub;
      }

      out.attr("class") = "mesh3d";
      return out;

    } catch (std::exception& e) {
      Rcpp::stop( e.what());
    } catch (...) {
      Rcpp::stop("unknown exception");
    }
  };

  static std::vector<bool> checkListNames(Rcpp::List mylist, Rcpp::CharacterVector mychar) {
    try {
      Rcpp::CharacterVector nam = mylist.names();
      Rcpp::IntegerVector ind(Rf_match(nam, mychar, 0));
      Rcpp::LogicalVector log(ind);
      std::vector<bool> out = Rcpp::as<std::vector<bool> >(log);
      return out;
    } catch (std::exception& e) {
      Rcpp::stop( e.what());
    } catch (...) {
      Rcpp::stop("unknown exception");
    }
  }

  static void mesh3d2vcg(MeshType &m, SEXP mesh_, bool zerobegin = false,
                         bool readnormals = true,bool readfaces = true) {
    Rcpp::List mesh(mesh_);
    Rcpp::CharacterVector mychar = Rcpp::CharacterVector::create("vb","it","normals");
    std::vector<bool> test = checkListNames(mesh, mychar);

    for (int i = 0; i < 3; i++) {
      if (!test[i]) {
        std::string tmp = Rcpp::as<std::string>(mychar[i]);
        mesh[tmp] = Rcpp::wrap(0);
      }
    }
    if (!test[0]) {
      Rcpp::stop("mesh has no vertices");
    }
    int out = vcgReadR(m, mesh["vb"], mesh["it"], mesh["normals"], zerobegin, readnormals, readfaces);

    if( out != 0 ) {
      Rcpp::stop("Unable to convert mesh");
    }
  };

};


// Additional definitions
typedef typename vcg::tri::UpdateTopology<MyMesh>::PEdge MyPEdge;
// Map from edge to new‐vertex
typedef std::pair<VertexPointer,VertexPointer> EdgePair;

struct EdgeHash {
  size_t operator()(EdgePair const &p) const noexcept
  {
    return std::hash<void*>()((void*)p.first) ^ (std::hash<void*>()((void*)p.second) << 1);
  }
};

struct EdgeEqual {
  bool operator()(
      EdgePair const &a,
      EdgePair const &b) const noexcept
  {
    return a.first==b.first && a.second==b.second;
  }
};

}

#endif // RAVETOOLS_VCG_COMMON_H
