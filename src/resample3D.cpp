#include <cmath>
#include "utils.h"
#include "TinyParallel.h"

using namespace Rcpp;


/*** R
# DIPSAUS DEBUG START
devtools::load_all()
arrayDim = rev(c(256, 128, 16))
fromArray = array(rnorm(prod(arrayDim)), rev(arrayDim))
oldVoxToWorld <- matrix(nrow = 4, byrow = TRUE, c(
  1.5, 0, 0.8, -200,
  0.01, 0, -1, 10,
  0, -0.8, 0.1, 60,
  0, 0, 0, 1
))
newVoxToWorld <- matrix(nrow = 4, byrow = TRUE, c(
  15, 0, 0.08, -140,
  0.1, -0.5, 0, 30,
  0, 0.1, 0.1, -13.6,
  0, 0, 0, 1
))

solve(oldVoxToWorld) %*% newVoxToWorld %*% c(arrayDim / 2,1)

reList <- resample3D(arrayDim, fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = 0.0)

dim(reList[[2]]) <- c(4,4)

mat <- solve(reList[[2]]) %*% solve(oldVoxToWorld) %*% newVoxToWorld
diag(mat) <- diag(mat) - 1
max(abs(mat))


# validate
sum(reList[[1]])

expected_array <- array(0, arrayDim)
idx <- t(cbind(arrayInd(seq_len(prod(arrayDim)), arrayDim) - 1, 1))
idx <- round(solve(oldVoxToWorld) %*% newVoxToWorld %*% idx)[1:3, , drop = FALSE]
idx[idx < 0 | idx >= dim(fromArray)] <- NA
idx <- colSums(idx * c(1, cumprod(dim(fromArray)))[1:3]) + 1
expected_array[!is.na(idx)] <- fromArray[idx[!is.na(idx)]]

range(reList[[1]] - expected_array)

# profile
arrayDim = rev(c(256, 128, 16))
fromArray = array(rnorm(prod(arrayDim)), rev(arrayDim))
oldVoxToWorld <- matrix(nrow = 4, byrow = TRUE, c(
  1.5, 0, 0.8, -200,
  0.01, 0, -1, 10,
  0, -0.8, 0.1, 60,
  0, 0, 0, 1
))
newVoxToWorld <- matrix(nrow = 4, byrow = TRUE, c(
  15, 0, 0.08, -140,
  0.1, -0.5, 0, 30,
  0, 0.1, 0.1, -13.6,
  0, 0, 0, 1
))

microbenchmark::microbenchmark(
  {
    ravetools::ravetools_threads(1)
    reList <- resample3D(rep(512, 3), fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = 0)
  },
  {
    ravetools::ravetools_threads(4)
    reList <- resample3D(rep(512, 3), fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = 0)
  },
  {
    ravetools::ravetools_threads(8)
    reList <- resample3D(rep(512, 3), fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = 0)
  }, times = 3, setup = { gc() }
)

*/


template <typename T>
struct Resampler3D : public TinyParallel::Worker
{
  // dimension of new array
  const R_xlen_t& nd1, nd2, nd3;

  // dimension of old array
  const R_xlen_t& od1, od2, od3;

  // vox from new to old, 3x4 matrix
  const double& a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34;

  const T& na_fill;

  T* &re_ptr;
  T* &x_ptr;

  // temporary, used to store cumprod of c(nd1, 2, 3)
  R_xlen_t cnd2, cnd3;

  Resampler3D(
    T* &x_ptr,
    T* &re_ptr,
    const T& na_fill,
    const R_xlen_t& nd1, const R_xlen_t& nd2, const R_xlen_t& nd3,
    const R_xlen_t& od1, const R_xlen_t& od2, const R_xlen_t& od3,
    const double& a11, const double& a12, const double& a13, const double& a14,
    const double& a21, const double& a22, const double& a23, const double& a24,
    const double& a31, const double& a32, const double& a33, const double& a34
  ):  nd1(nd1), nd2(nd2), nd3(nd3), od1(od1), od2(od2), od3(od3),
      a11 (a11), a12 (a12), a13 (a13), a14 (a14),
      a21 (a21), a22 (a22), a23 (a23), a24 (a24),
      a31 (a31), a32 (a32), a33 (a33), a34 (a34),
      na_fill(na_fill), re_ptr (re_ptr), x_ptr (x_ptr)
  {
    this->cnd2 = nd1;
    this->cnd3 = nd1 * nd2;
  }

  void operator()(std::size_t begin, std::size_t end) {
    R_xlen_t begin_ = (R_xlen_t) begin;
    R_xlen_t end_ = (R_xlen_t) end;

    // for transforming indices
    R_xlen_t vi_int, vj_int, vk_int;
    double vi_f, vj_f, vk_f;

    for( R_xlen_t ii = begin_; ii < end_; ii++ ) {
      // get IJK for the new array
      vk_int = ii / this->cnd3;
      vi_int = ii - vk_int * this->cnd3;
      vj_int = vi_int / this->cnd2;
      vi_int -= vj_int * this->cnd2;

      // Rcout << ii << " " << vi_int << " " << vj_int << " " << vk_int << "-----\n";

      // get voxel index in old array
      vi_f = std::nearbyint( this->a11 * vi_int + this->a12 * vj_int + this->a13 * vk_int + this->a14 );
      vj_f = std::nearbyint( this->a21 * vi_int + this->a22 * vj_int + this->a23 * vk_int + this->a24 );
      vk_f = std::nearbyint( this->a31 * vi_int + this->a32 * vj_int + this->a33 * vk_int + this->a34 );

      vi_int = (R_xlen_t) vi_f;
      vj_int = (R_xlen_t) vj_f;
      vk_int = (R_xlen_t) vk_f;

      // Rcout << ii << " " << vi_int << " " << vj_int << " " << vk_int << " " << (vi_int + od1 * ( vj_int + od2 * vk_int )) << "\n";

      if( vi_int < 0 || vi_int >= this->od1 || vj_int < 0 || vj_int >= this->od2 || vk_int < 0 || vk_int >= this->od3 ) {
        // outofbound, use NA values
        *(this->re_ptr + ii) = this->na_fill;
      } else {
        // Rcout << ii << " " << vi_int << "\n";
        *(this->re_ptr + ii) = *(this->x_ptr + (vi_int + this->od1 * ( vj_int + this->od2 * vk_int )));
      }
      // if( vi_int == 0 && vj_int == 0 ) {
      //   Rcpp::checkUserInterrupt();
      // }
    }

  }
};




/**
 * Resample a `arrayDim` array from `fromArray`, with new
 * vox2ras=newVoxToWorld and old vox2ras=oldVoxToWorld.
 * The resampling happens in the world space.
 *
 * newVoxToWorld and oldVoxToWorld should have at least 12 elements
 * since R is col-major, the easiest way to avoid checking the last 4 elems
 * is to transpose the matrix and only examine the first 12 elements
 * hence newVoxToWorldTransposed = t(newVoxToWorld), and
 * oldVoxToWorldTransposed = t(oldVoxToWorld)
 */
// [[Rcpp::export]]
SEXP resample3D(const SEXP& arrayDim,
                const SEXP& fromArray,
                const SEXP& newVoxToWorldTransposed,
                const SEXP& oldVoxToWorldTransposed,
                const SEXP& na) {

  if( XLENGTH(arrayDim) < 3 ) {
    Rcpp::stop("C++ `resample3D`: the dimension of new array must be 3D.");
  }
  SEXP oldDim = Rf_getAttrib(fromArray, R_DimSymbol);

  if( XLENGTH(oldDim) < 3 ) {
    Rcpp::stop("C++ `resample3D`: the dimension of sampling volume must be 3D.");
  }

  // gather dimension information
  R_xlen_t od1, od2, od3, nd1, nd2, nd3;
  if(TYPEOF(oldDim) == INTSXP) {
    int* ptr = INTEGER(oldDim);
    od1 = ptr[0];
    od2 = ptr[1];
    od3 = ptr[2];
  } else {
    SEXP oldDim_ = PROTECT(Rf_coerceVector(oldDim, INTSXP));
    int* ptr = INTEGER(oldDim_);
    od1 = ptr[0];
    od2 = ptr[1];
    od3 = ptr[2];
    UNPROTECT(1); // oldDim_
  }

  if(TYPEOF(arrayDim) == INTSXP) {
    int* ptr = INTEGER(arrayDim);
    nd1 = ptr[0];
    nd2 = ptr[1];
    nd3 = ptr[2];
  } else {
    SEXP arrayDim_ = PROTECT(Rf_coerceVector(arrayDim, INTSXP));
    int* ptr = INTEGER(arrayDim_);
    nd1 = ptr[0];
    nd2 = ptr[1];
    nd3 = ptr[2];
    UNPROTECT(1); // arrayDim_
  }


  // check if newVoxToWorldTransposed and newVoxToWorldTransposed are matrices
  if(XLENGTH(newVoxToWorldTransposed) < 12) {
    Rcpp::stop("C++ `resample3D`: the voxel-to-world matrix for the sampling volume is invalid.");
  }
  if(XLENGTH(oldVoxToWorldTransposed) < 12) {
    Rcpp::stop("C++ `resample3D`: the voxel-to-world matrix for new array is invalid.");
  }

  double  b11, b12, b13, b14,
          b21, b22, b23, b24,
          b31, b32, b33, b34,
          b41 = 0.0, b42 = 0.0, b43 = 0.0, b44 = 1.0;

  if( TYPEOF(newVoxToWorldTransposed) != REALSXP ) {
    SEXP vox2rasNew = PROTECT(Rf_coerceVector(newVoxToWorldTransposed, REALSXP));
    double* ptr = REAL(vox2rasNew);
    b11 = *ptr++; b12 = *ptr++; b13 = *ptr++; b14 = *ptr++;
    b21 = *ptr++; b22 = *ptr++; b23 = *ptr++; b24 = *ptr++;
    b31 = *ptr++; b32 = *ptr++; b33 = *ptr++; b34 = *ptr++;
    UNPROTECT(1);
  } else {
    double* ptr = REAL(newVoxToWorldTransposed);
    b11 = *ptr++; b12 = *ptr++; b13 = *ptr++; b14 = *ptr++;
    b21 = *ptr++; b22 = *ptr++; b23 = *ptr++; b24 = *ptr++;
    b31 = *ptr++; b32 = *ptr++; b33 = *ptr++; b34 = *ptr++;
  }

  double  a11, a12, a13, a14,
          a21, a22, a23, a24,
          a31, a32, a33, a34,
          a41 = 0.0, a42 = 0.0, a43 = 0.0, a44 = 1.0;

  if( TYPEOF(oldVoxToWorldTransposed) != REALSXP ) {
    SEXP vox2rasOld = PROTECT(Rf_coerceVector(oldVoxToWorldTransposed, REALSXP));
    double* ptr = REAL(vox2rasOld);
    a11 = *ptr++; a12 = *ptr++; a13 = *ptr++; a14 = *ptr++;
    a21 = *ptr++; a22 = *ptr++; a23 = *ptr++; a24 = *ptr++;
    a31 = *ptr++; a32 = *ptr++; a33 = *ptr++; a34 = *ptr++;
    UNPROTECT(1);
  } else {
    double* ptr = REAL(oldVoxToWorldTransposed);
    a11 = *ptr++; a12 = *ptr++; a13 = *ptr++; a14 = *ptr++;
    a21 = *ptr++; a22 = *ptr++; a23 = *ptr++; a24 = *ptr++;
    a31 = *ptr++; a32 = *ptr++; a33 = *ptr++; a34 = *ptr++;
  }

  // calculate transform from new vox to old vox
  // a is old vox2ras, b is new
  // solve(a) %*% b is from new vox to old vox
  double  t11 = a23 * a34 * a42 - a24 * a33 * a42 + a24 * a32 * a43 - a22 * a34 * a43 - a23 * a32 * a44 + a22 * a33 * a44,
          t12 = a14 * a33 * a42 - a13 * a34 * a42 - a14 * a32 * a43 + a12 * a34 * a43 + a13 * a32 * a44 - a12 * a33 * a44,
          t13 = a13 * a24 * a42 - a14 * a23 * a42 + a14 * a22 * a43 - a12 * a24 * a43 - a13 * a22 * a44 + a12 * a23 * a44,
          t14 = a14 * a23 * a32 - a13 * a24 * a32 - a14 * a22 * a33 + a12 * a24 * a33 + a13 * a22 * a34 - a12 * a23 * a34;

  double det = a11 * t11 + a21 * t12 + a31 * t13 + a41 * t14;

  if ( det < 0.0000001 && det > -0.0000001 ) {
    a11 = a12 = a13 = a14 = a21 = a22 = a23 = a24 = a31 = a32 = a33 = a34 = 0.0;
  } else {
    double detInv = 1.0 / det;
    double  n11 = t11 * detInv,
            n21 = ( a24 * a33 * a41 - a23 * a34 * a41 - a24 * a31 * a43 + a21 * a34 * a43 + a23 * a31 * a44 - a21 * a33 * a44 ) * detInv,
            n31 = ( a22 * a34 * a41 - a24 * a32 * a41 + a24 * a31 * a42 - a21 * a34 * a42 - a22 * a31 * a44 + a21 * a32 * a44 ) * detInv,
            n41 = ( a23 * a32 * a41 - a22 * a33 * a41 - a23 * a31 * a42 + a21 * a33 * a42 + a22 * a31 * a43 - a21 * a32 * a43 ) * detInv,

            n12 = t12 * detInv,
            n22 = ( a13 * a34 * a41 - a14 * a33 * a41 + a14 * a31 * a43 - a11 * a34 * a43 - a13 * a31 * a44 + a11 * a33 * a44 ) * detInv,
            n32 = ( a14 * a32 * a41 - a12 * a34 * a41 - a14 * a31 * a42 + a11 * a34 * a42 + a12 * a31 * a44 - a11 * a32 * a44 ) * detInv,
            n42 = ( a12 * a33 * a41 - a13 * a32 * a41 + a13 * a31 * a42 - a11 * a33 * a42 - a12 * a31 * a43 + a11 * a32 * a43 ) * detInv,

            n13 = t13 * detInv,
            n23 = ( a14 * a23 * a41 - a13 * a24 * a41 - a14 * a21 * a43 + a11 * a24 * a43 + a13 * a21 * a44 - a11 * a23 * a44 ) * detInv,
            n33 = ( a12 * a24 * a41 - a14 * a22 * a41 + a14 * a21 * a42 - a11 * a24 * a42 - a12 * a21 * a44 + a11 * a22 * a44 ) * detInv,
            n43 = ( a13 * a22 * a41 - a12 * a23 * a41 - a13 * a21 * a42 + a11 * a23 * a42 + a12 * a21 * a43 - a11 * a22 * a43 ) * detInv,

            n14 = t14 * detInv,
            n24 = ( a13 * a24 * a31 - a14 * a23 * a31 + a14 * a21 * a33 - a11 * a24 * a33 - a13 * a21 * a34 + a11 * a23 * a34 ) * detInv,
            n34 = ( a14 * a22 * a31 - a12 * a24 * a31 - a14 * a21 * a32 + a11 * a24 * a32 + a12 * a21 * a34 - a11 * a22 * a34 ) * detInv,
            n44 = ( a12 * a23 * a31 - a13 * a22 * a31 + a13 * a21 * a32 - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33 ) * detInv;

    a11 = n11 * b11 + n12 * b21 + n13 * b31 + n14 * b41;
    a12 = n11 * b12 + n12 * b22 + n13 * b32 + n14 * b42;
    a13 = n11 * b13 + n12 * b23 + n13 * b33 + n14 * b43;
    a14 = n11 * b14 + n12 * b24 + n13 * b34 + n14 * b44;

    a21 = n21 * b11 + n22 * b21 + n23 * b31 + n24 * b41;
    a22 = n21 * b12 + n22 * b22 + n23 * b32 + n24 * b42;
    a23 = n21 * b13 + n22 * b23 + n23 * b33 + n24 * b43;
    a24 = n21 * b14 + n22 * b24 + n23 * b34 + n24 * b44;

    a31 = n31 * b11 + n32 * b21 + n33 * b31 + n34 * b41;
    a32 = n31 * b12 + n32 * b22 + n33 * b32 + n34 * b42;
    a33 = n31 * b13 + n32 * b23 + n33 * b33 + n34 * b43;
    a34 = n31 * b14 + n32 * b24 + n33 * b34 + n34 * b44;

    a41 = n41 * b11 + n42 * b21 + n43 * b31 + n44 * b41;
    a42 = n41 * b12 + n42 * b22 + n43 * b32 + n44 * b42;
    a43 = n41 * b13 + n42 * b23 + n43 * b33 + n44 * b43;
    a44 = n41 * b14 + n42 * b24 + n43 * b34 + n44 * b44;

  }

  // now matrix a is the transform from new vox to old vox

  // construct new array
  R_xlen_t retLen = nd1 * nd2 * nd3;

  unsigned int arrayType = (unsigned int) TYPEOF(fromArray);
  SEXP re = PROTECT(Rf_allocVector(arrayType, retLen));

  // double* re_ptr = REAL(re);
  // double* x_ptr = REAL(oldArray);
  // R_xlen_t cnd2 = nd1, cnd3 = nd1 * nd2;
  // R_xlen_t vi_int, vj_int, vk_int;
  // float vi_f, vj_f, vk_f;
  // for( R_xlen_t ii = 0; ii < retLen; ii++ ) {
  //   // get IJK for the new array
  //   vk_int = ii / cnd3;
  //   vi_int = ii - vk_int * cnd3;
  //   vj_int = vi_int / cnd2;
  //   vi_int -= vj_int * cnd2;
  //
  //   // Rcout << ii << " " << vi_int << " " << vj_int << " " << vk_int << "-----\n";
  //
  //   // get voxel index in old array
  //   vi_f = std::roundf( a11 * vi_int + a12 * vj_int + a13 * vk_int + a14 );
  //   vj_f = std::roundf( a21 * vi_int + a22 * vj_int + a23 * vk_int + a24 );
  //   vk_f = std::roundf( a31 * vi_int + a32 * vj_int + a33 * vk_int + a34 );
  //
  //   vi_int = (R_xlen_t) vi_f;
  //   vj_int = (R_xlen_t) vj_f;
  //   vk_int = (R_xlen_t) vk_f;
  //
  //   // Rcout << ii << " " << vi_int << " " << vj_int << " " << vk_int << " " << (vi_int + od1 * ( vj_int + od2 * vk_int )) << "\n";
  //
  //   if( vi_int < 0 || vi_int >= od1 || vj_int < 0 || vj_int >= od2 || vk_int < 0 || vk_int >= od3 ) {
  //     // outofbound, use NA values
  //     *(re_ptr + ii) = na;
  //   } else {
  //     // Rcout << ii << " " << vi_int << "\n";
  //     *(re_ptr + ii) = *(x_ptr + (vi_int + od1 * ( vj_int + od2 * vk_int )));
  //   }
  //   // if( vi_int == 0 && vj_int == 0 ) {
  //   //   Rcpp::checkUserInterrupt();
  //   // }
  // }

  switch(arrayType) {

    case REALSXP: {
      const double na_ = REAL(na)[0];
      double* x_ptr = REAL(fromArray);
      double* re_ptr = REAL(re);
      Resampler3D<double> sampler(
          x_ptr, re_ptr, na_,
          nd1, nd2, nd3, od1, od2, od3,
          a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34);
      parallelFor(0, retLen, sampler);
      break;
    }
    case INTSXP: {
      const int na_ = INTEGER(na)[0];
      int* x_ptr = INTEGER(fromArray);
      int* re_ptr = INTEGER(re);
      Resampler3D<int> sampler(
          x_ptr, re_ptr, na_,
          nd1, nd2, nd3, od1, od2, od3,
          a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34);
      parallelFor(0, retLen, sampler);
      break;
    }
    case LGLSXP: {
      const int na_ = LOGICAL(na)[0];
      int* x_ptr = LOGICAL(fromArray);
      int* re_ptr = LOGICAL(re);
      Resampler3D<int> sampler(
          x_ptr, re_ptr, na_,
          nd1, nd2, nd3, od1, od2, od3,
          a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34);
      parallelFor(0, retLen, sampler);
      break;
    }
    case CPLXSXP: {
      const Rcomplex na_ = COMPLEX(na)[0];
      Rcomplex* x_ptr = COMPLEX(fromArray);
      Rcomplex* re_ptr = COMPLEX(re);
      Resampler3D<Rcomplex> sampler(
          x_ptr, re_ptr, na_,
          nd1, nd2, nd3, od1, od2, od3,
          a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34);
      parallelFor(0, retLen, sampler);
      break;
    }

    case RAWSXP: {
      const Rbyte na_ = RAW(na)[0];
      Rbyte* x_ptr = RAW(fromArray);
      Rbyte* re_ptr = RAW(re);
      Resampler3D<Rbyte> sampler(
          x_ptr, re_ptr, na_,
          nd1, nd2, nd3, od1, od2, od3,
          a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34);
      parallelFor(0, retLen, sampler);
      break;
    }

    case STRSXP: {
      // int64 is the same size as double, with 8 bytes
      const int64_t na_ = -1;

      SEXP re_ = PROTECT(Rf_allocVector(REALSXP, retLen));
      int64_t* re_ptr = (int64_t*) REAL(re_);

      int64_t odLen = (int64_t) (od1 * od2 * od3);
      SEXP fromArray_ = PROTECT(Rf_allocVector(REALSXP, odLen));

      int64_t* x_ptr = (int64_t*) REAL(fromArray_);
      for(int64_t ii = 0 ; ii < odLen; ii++) {
        *x_ptr++ = ii;
      }
      x_ptr = (int64_t*) REAL(fromArray_);

      Resampler3D<int64_t> sampler(
          x_ptr, re_ptr, na_,
          nd1, nd2, nd3, od1, od2, od3,
          a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34);
      parallelFor(0, retLen, sampler);

      UNPROTECT(1); // fromArray_

      const SEXP na_str = STRING_ELT(na, 0);
      re_ptr = (int64_t*) REAL(re_);
      for(int64_t ii = 0 ; ii < retLen; ii++) {
        const int64_t& reIdx = *re_ptr++;
        if( reIdx == na_ ) {
          SET_STRING_ELT(re, ii, na_str);
        } else {
          const SEXP& el = STRING_ELT(fromArray, reIdx);
          SET_STRING_ELT(re, ii, el);
        }
      }

      UNPROTECT(1); // re_
      break;
    }

    default: {
      Rcpp::stop("C++ `resample3D`: Unable to resample volume due to unrecognized storage format.");
    }

  }

  // set dimensions
  SEXP newDim_ = PROTECT(Rf_allocVector(INTSXP, 3));

  int* newDim_ptr = INTEGER(newDim_);
  *newDim_ptr = (int) nd1;
  *(newDim_ptr + 1) = (int) nd2;
  *(newDim_ptr + 2) = (int) nd3;

  Rf_setAttrib(re, R_DimSymbol, newDim_);

  // save transforms
  SEXP vox2vox = PROTECT(Rf_allocVector(REALSXP, 16));
  double* vox2vox_ptr = REAL(vox2vox);
  vox2vox_ptr[0] = a11;
  vox2vox_ptr[1] = a21;
  vox2vox_ptr[2] = a31;
  vox2vox_ptr[3] = a41;

  vox2vox_ptr[4] = a12;
  vox2vox_ptr[5] = a22;
  vox2vox_ptr[6] = a32;
  vox2vox_ptr[7] = a42;

  vox2vox_ptr[8] = a13;
  vox2vox_ptr[9] = a23;
  vox2vox_ptr[10] = a33;
  vox2vox_ptr[11] = a43;

  vox2vox_ptr[12] = a14;
  vox2vox_ptr[13] = a24;
  vox2vox_ptr[14] = a34;
  vox2vox_ptr[15] = a44;

  SEXP reList = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(reList, 0, re);
  SET_VECTOR_ELT(reList, 1, vox2vox);


  UNPROTECT(4); // reList, vox2vox, newDim_, re

  return(reList);
}



