#include <Rcpp.h>
#include "vcgCommon.h"

// [[Rcpp::export]]
SEXP vcgIsoSurface(SEXP array_, double thresh) {
  try {
    Rcpp::IntegerVector arrayDims( Rf_getAttrib(array_, R_DimSymbol) );
    std::vector<float> vecArray = Rcpp::as<std::vector<float> >(array_);

    ravetools::MyMesh m;
    ravetools::VertexIterator vi;
    ravetools::FaceIterator fi;
    int i,j,k;

    // typedef MySimpleVolume<ravetools::MySimpleVoxel> MyVolume;
    ravetools::MyVolume	volume;
    // typedef vcg::tri::TrivialWalker<ravetools::MyMesh, ravetools::MyVolume>	MyWalker;
    // typedef vcg::tri::MarchingCubes<ravetools::MyMesh, ravetools::MyWalker>	MyMarchingCubes;
    ravetools::MyWalker walker;
    volume.Init( vcg::Point3i(arrayDims[0], arrayDims[1], arrayDims[2]) );
    for( i = 0; i < arrayDims[0]; i++ ) {
      for( j = 0; j < arrayDims[1]; j++ ) {
        for( k = 0; k < arrayDims[2]; k++ ) {
          int tmpval = vecArray[ i + j * arrayDims[0] + k * ( arrayDims[0] * arrayDims[1] ) ];
          /*if (tmpval >= lower && tmpval <= upper)
           volume.Val(i,j,k)=tmpval;
           else*/
          volume.Val(i,j,k)=tmpval;
        }
      }
    }

    Rcpp::checkUserInterrupt();
    //write back
    /*volume.Init(Point3i(64,64,64));
     for(int i=0;i<64;i++)
     for(int j=0;j<64;j++)
     for(int k=0;k<64;k++)
     volume.Val(i,j,k)=(j-32)*(j-32)+(k-32)*(k-32)  + i*10*(float)math::Perlin::Noise(i*.2,j*.2,k*.2);*/
    ravetools::MyMarchingCubes	mc(m, walker);
    walker.BuildMesh<ravetools::MyMarchingCubes>(m, volume, mc, thresh);
    vcg::tri::Allocator< ravetools::MyMesh >::CompactVertexVector(m);
    vcg::tri::Allocator< ravetools::MyMesh >::CompactFaceVector(m);
    vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal< ravetools::MyMesh >::NormalizePerVertex(m);
    vcg::SimpleTempData< ravetools::MyMesh::VertContainer, int > indiceout(m.vert);
    Rcpp::NumericMatrix vbout(3,m.vn), normals(3,m.vn);
    Rcpp::IntegerMatrix itout(3,m.fn);

    Rcpp::checkUserInterrupt();

    vi=m.vert.begin();
    for (i=0;  i < m.vn; i++) {
      indiceout[vi] = i;
      vbout(0,i) = (*vi).P()[0];
      vbout(1,i) = (*vi).P()[1];
      vbout(2,i) = (*vi).P()[2];
      normals(0,i) = (*vi).N()[0];
      normals(1,i) = (*vi).N()[1];
      normals(2,i) = (*vi).N()[2];
      ++vi;
    }
    ravetools::FacePointer fp;

    fi=m.face.begin();
    j = 0;
    for (i=0; i < m.fn; i++) {
      fp=&(*fi);
      itout(0,i) = indiceout[fp->cV(0)]+1;
      itout(1,i) = indiceout[fp->cV(1)]+1;
      itout(2,i) = indiceout[fp->cV(2)]+1;
      ++fi;
    }
    //delete &walker;
    //delete &volume;
    return Rcpp::List::create(Rcpp::Named("vb") = vbout,
                              Rcpp::Named("it") = itout,
                              Rcpp::Named("normals") = normals);
  } catch (std::exception& e) {
    Rcpp::stop( e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgSmoothImplicit(
    SEXP vb_, SEXP it_, float lambda_, bool useMassMatrix, bool fixBorder,
    bool useCotWeight, int degree, float lapWeight_, bool SmoothQ)
{
  try {
    int i;
    ravetools::MyMesh m;
    ravetools::VertexIterator vi;
    ravetools::FaceIterator fi;

    ravetools::ScalarType lambda = lambda_;
    ravetools::ScalarType lapWeight = lapWeight_;

    //allocate mesh and fill it
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m,vb_,it_);


    vcg::ImplicitSmoother<ravetools::MyMesh>::Parameter par;
    par.lambda = lambda;
    par.useMassMatrix = useMassMatrix;
    par.fixBorder = fixBorder;
    par.useCotWeight = useCotWeight;
    par.degree = degree;
    par.lapWeight = lapWeight;
    par.SmoothQ = SmoothQ;

    vcg::ImplicitSmoother<ravetools::MyMesh>::Compute(m, par);

    Rcpp::checkUserInterrupt();

    vcg::tri::Allocator<ravetools::MyMesh>::CompactVertexVector(m);
    vcg::tri::Allocator<ravetools::MyMesh>::CompactFaceVector(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(m);
    Rcpp::NumericMatrix vb(3, m.vn);
    Rcpp::NumericMatrix normals(3, m.vn);
    Rcpp::IntegerMatrix itout(3, m.fn);
    //write back output
    vcg::SimpleTempData<ravetools::MyMesh::VertContainer,int>indices(m.vert);

    Rcpp::checkUserInterrupt();

    // write back updated mesh
    vi=m.vert.begin();
    for (i=0; i < m.vn; i++) {
      indices[vi] = i;
      if( ! vi->IsD() ) {
        vb(0,i) = (*vi).P()[0];
        vb(1,i) = (*vi).P()[1];
        vb(2,i) = (*vi).P()[2];
        normals(0,i) = (*vi).N()[0];
        normals(1,i) = (*vi).N()[1];
        normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }

    ravetools::FacePointer fp;
    fi=m.face.begin();
    for (i=0; i < m.fn; i++) {
      fp=&(*fi);
      if( ! fp->IsD() ) {
        itout(0,i) = indices[fp->cV(0)]+1;
        itout(1,i) = indices[fp->cV(1)]+1;
        itout(2,i) = indices[fp->cV(2)]+1;
      }
      ++fi;
    }
    return Rcpp::List::create(Rcpp::Named("vb") = vb,
                              Rcpp::Named("normals") = normals,
                              Rcpp::Named("it") = itout
    );

  } catch (std::exception& e) {
    Rcpp::stop( e.what());
    return Rcpp::wrap(1);
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}


// [[Rcpp::export]]
SEXP vcgSmooth(SEXP vb_, SEXP it_, int iter, int method, float lambda, float mu, float delta_)
{
  try {
    int i;
    ravetools::MyMesh m;
    ravetools::VertexIterator vi;
    ravetools::FaceIterator fi;
    //set up parameters
    ravetools::ScalarType delta = delta_;
    //allocate mesh and fill it
    ravetools::IOMesh<ravetools::MyMesh>::vcgReadR(m,vb_,it_);

    Rcpp::checkUserInterrupt();

    if (method == 0) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      unsigned int cnt = vcg::tri::UpdateSelection<ravetools::MyMesh>::VertexFromFaceStrict(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordTaubin(m, iter, lambda, mu, cnt>0);
    } else if (method == 1) {
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordLaplacian(m, iter);
    } else if (method == 2) {
      vcg::tri::UpdateSelection<ravetools::MyMesh>::FaceAll(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      unsigned int cnt=vcg::tri::UpdateSelection<ravetools::MyMesh>::VertexFromFaceStrict(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordLaplacianHC(m, iter,cnt>0);
    } else if (method == 3) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceClearB(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordScaleDependentLaplacian_Fujiwara(m,iter,delta);
    } else if (method == 4) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceClearB(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordLaplacianAngleWeighted(m,iter,delta);
    }
    else if (method == 5) {
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceBorderFromNone(m);
      vcg::tri::UpdateFlags<ravetools::MyMesh>::FaceClearB(m);
      vcg::tri::Smooth<ravetools::MyMesh>::VertexCoordPlanarLaplacian(m, iter, delta);
    }

    Rcpp::checkUserInterrupt();

    vcg::tri::Allocator<ravetools::MyMesh>::CompactVertexVector(m);
    vcg::tri::Allocator<ravetools::MyMesh>::CompactFaceVector(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexAngleWeighted(m);
    vcg::tri::UpdateNormal<ravetools::MyMesh>::NormalizePerVertex(m);
    Rcpp::NumericMatrix vb(3, m.vn);
    Rcpp::NumericMatrix normals(3, m.vn);
    Rcpp::IntegerMatrix itout(3, m.fn);
    //write back output
    vcg::SimpleTempData<ravetools::MyMesh::VertContainer,int>indices(m.vert);

    // write back updated mesh
    vi=m.vert.begin();
    for (i=0; i < m.vn; i++) {
      indices[vi] = i;
      if( ! vi->IsD() ) {
        vb(0,i) = (*vi).P()[0];
        vb(1,i) = (*vi).P()[1];
        vb(2,i) = (*vi).P()[2];
        normals(0,i) = (*vi).N()[0];
        normals(1,i) = (*vi).N()[1];
        normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }

    ravetools::FacePointer fp;
    fi=m.face.begin();
    for (i=0; i < m.fn; i++) {
      fp=&(*fi);
      if( ! fp->IsD() ) {
        itout(0,i) = indices[fp->cV(0)]+1;
        itout(1,i) = indices[fp->cV(1)]+1;
        itout(2,i) = indices[fp->cV(2)]+1;
      }
      ++fi;
    }
    return Rcpp::List::create(Rcpp::Named("vb") = vb,
                              Rcpp::Named("normals") = normals,
                              Rcpp::Named("it") = itout
    );

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
    return Rcpp::wrap(1);
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgUniformResample(
    const SEXP& vb_, const SEXP& it_, const float& voxelSize, const float& offsetThr,
    const bool& discretizeFlag, const bool& multiSampleFlag, const bool& absDistFlag,
    const bool& mergeCloseVert, const bool& silent) {
  try {
    ravetools::MyMesh m, baseMesh, offsetMesh;
    ravetools::IOMesh< ravetools::MyMesh >::vcgReadR( baseMesh , vb_ , it_ );
    if ( baseMesh.fn == 0 ) {
      Rcpp::stop( "This filter requires a mesh with some faces, it does not work on point cloud");
    }
    vcg::tri::UpdateBounding< ravetools::MyMesh >::Box(baseMesh);
    baseMesh.face.EnableNormal();
    vcg::Point3i volumeDim;
    vcg::Box3f volumeBox = baseMesh.bbox;
    volumeBox.Offset( volumeBox.Diag()/10.0f + offsetThr );

    BestDim(volumeBox , voxelSize, volumeDim );

    Rcpp::checkUserInterrupt();

    if (!silent) {
      Rprintf("Resampling mesh using a volume of %i x %i x %i\n", volumeDim[0], volumeDim[1], volumeDim[2] );
      Rprintf("  VoxelSize is %f, offset is %f\n", voxelSize, offsetThr);
      Rprintf("  Mesh Box is %f %f %f\n", baseMesh.bbox.DimX(),
              baseMesh.bbox.DimY(), baseMesh.bbox.DimZ() );
    }
    vcg::tri::Resampler< ravetools::MyMesh, ravetools::MyMesh >::Resample(
        baseMesh, offsetMesh, volumeBox, volumeDim, voxelSize * 3.5f,
        offsetThr, discretizeFlag, multiSampleFlag, absDistFlag
    );
    Rcpp::checkUserInterrupt();
    if ( mergeCloseVert ) {
      float mergeThr = offsetMesh.bbox.Diag() / 10000.0f;
      int total = vcg::tri::Clean< ravetools::MyMesh >::MergeCloseVertex( offsetMesh , mergeThr );
      if ( !silent ) {
        Rprintf("Successfully merged %d vertices with a distance lower than %f\n", total, mergeThr);
      }
    }
    vcg::tri::Allocator< ravetools::MyMesh >::CompactVertexVector( offsetMesh );
    vcg::tri::Allocator< ravetools::MyMesh >::CompactFaceVector( offsetMesh );
    vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertexAngleWeighted( offsetMesh );
    vcg::tri::UpdateNormal< ravetools::MyMesh >::NormalizePerVertex( offsetMesh );
    Rcpp::NumericMatrix vbout(3, offsetMesh.vn), normals(3, offsetMesh.vn);
    Rcpp::IntegerMatrix itout(3, offsetMesh.fn);
    vcg::SimpleTempData< ravetools::MyMesh::VertContainer, int > indiceout( offsetMesh.vert );
    ravetools::VertexIterator vi;
    vi = offsetMesh.vert.begin();
    for ( int i = 0 ; i < offsetMesh.vn ; i++ ) {
      indiceout[vi] = i;
      vbout(0,i) = (*vi).P()[0];
      vbout(1,i) = (*vi).P()[1];
      vbout(2,i) = (*vi).P()[2];
      normals(0,i) = (*vi).N()[0];
      normals(1,i) = (*vi).N()[1];
      normals(2,i) = (*vi).N()[2];
      ++vi;
    }
    ravetools::FaceIterator fi = offsetMesh.face.begin();
    for ( int i = 0; i < offsetMesh.fn ; i++, fi++ ) {
      itout(0, i) = indiceout[ fi->cV(0) ] + 1;
      itout(1, i) = indiceout[ fi->cV(1) ] + 1;
      itout(2, i) = indiceout[ fi->cV(2) ] + 1;
    }

    return Rcpp::List::create(Rcpp::Named("vb") = vbout,
                              Rcpp::Named("it") = itout,
                              Rcpp::Named("normals")=normals);

    return Rcpp::wrap(0);
  } catch ( std::exception& e ) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgUpdateNormals(SEXP vb_, SEXP it_, const int & select,
                      const Rcpp::IntegerVector & pointcloud, const bool & silent)
{
  try {
    ravetools::MyMesh m;
    ravetools::VertexIterator vi;

    // allocate mesh and fill it
    int check = ravetools::IOMesh< ravetools::MyMesh >::vcgReadR(m,vb_,it_);
    Rcpp::NumericMatrix normals(3, m.vn);
    if (check < 0) {
      Rcpp::stop("mesh has no faces and/or no vertices");
    } else if (check == 1) {
      if ( !silent ) {
        Rprintf("%s\n", "Info: mesh has no faces normals for point clouds are computed");
      }
      vcg::tri::PointCloudNormal< ravetools::MyMesh >::Param p;
      p.fittingAdjNum = pointcloud[0];
      p.smoothingIterNum = pointcloud[1];
      p.viewPoint = vcg::Point3f(0,0,0);
      p.useViewPoint = false;
      vcg::tri::PointCloudNormal< ravetools::MyMesh >::Compute(m,p);
    }  else {
      // update normals
      if (select == 0) {
        vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertex(m);
      } else {
        vcg::tri::UpdateNormal< ravetools::MyMesh >::PerVertexAngleWeighted(m);
      }
      vcg::tri::UpdateNormal< ravetools::MyMesh >::NormalizePerVertex(m);

      //write back
    }
    vi = m.vert.begin();
    vcg::SimpleTempData< ravetools::MyMesh::VertContainer , int > indiceout(m.vert);
    for ( int i = 0 ; i < m.vn ; i++) {
      if( ! vi->IsD() )	{
        normals(0,i) = (*vi).N()[0];
        normals(1,i) = (*vi).N()[1];
        normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }

    return Rcpp::wrap(normals);

  } catch (std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}



// [[Rcpp::export]]
SEXP vcgVolume( SEXP mesh_ )
{
  try {
    ravetools::MyMesh m;
    ravetools::IOMesh<ravetools::MyMesh>::mesh3d2vcg(m, mesh_);
    bool Watertight, Oriented = false;
    int VManifold, FManifold;
    float Volume = 0;
    // int numholes, BEdges = 0;
    //check manifoldness
    m.vert.EnableVFAdjacency();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    m.face.EnableNormal();
    vcg::tri::UpdateTopology<ravetools::MyMesh>::FaceFace(m);
    VManifold = vcg::tri::Clean<ravetools::MyMesh>::CountNonManifoldVertexFF(m);
    FManifold = vcg::tri::Clean<ravetools::MyMesh>::CountNonManifoldEdgeFF(m);

    if ((VManifold>0) || (FManifold>0)) {
      throw std::runtime_error(
        (
            "Mesh is not manifold\n  Non-manifold vertices: " +
              std::to_string(VManifold) +"\n" +
              "  Non-manifold edges: " +
              std::to_string(FManifold) +"\n"
        ).c_str()
      );
    }


    Watertight = vcg::tri::Clean<ravetools::MyMesh>::IsWaterTight(m);
    Oriented = vcg::tri::Clean<ravetools::MyMesh>::IsCoherentlyOrientedMesh(m);
    vcg::tri::Inertia<ravetools::MyMesh> mm(m);
    mm.Compute(m);
    Volume = mm.Mass();

    // the sign of the volume depend on the mesh orientation
    if (Volume < 0.0)
      Volume = -Volume;
    if (!Watertight)
      ::Rf_warning("Mesh is not watertight! USE RESULT WITH CARE!\n");
    if (!Oriented)
      ::Rf_warning("Mesh is not coherently oriented! USE RESULT WITH CARE!\n");

    return Rcpp::wrap(Volume);

  } catch (std::exception& e) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}


// [[Rcpp::export]]
SEXP vcgSphere(const int& subdiv, bool normals) {
  try {
    ravetools::MyMesh m;
    Sphere(m,subdiv);
    if (normals)
      vcg::tri::UpdateNormal<ravetools::MyMesh>::PerVertexNormalized(m);
    Rcpp::List out = ravetools::IOMesh<ravetools::MyMesh>::vcgToR(m,normals);
    return out;
  } catch (std::exception& e) {
    Rcpp::stop( e.what() );
  } catch (...) {
    Rcpp::stop("unknown exception");
  }
  return R_NilValue; // -Wall
}

