#include "fftw-wrapper.h"
#include "ffts.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fftw_r2c(SEXP data, int HermConj = 1,
              int fftwplanopt = 0,
              SEXP ret = R_NilValue,
              bool inplace = false) {
  int nprot = 0;

  // check HermConj and ret
  int xlen = Rf_length(data);
  int retlen = 0;
  if( HermConj == 1 ){
    retlen = xlen;
  } else {
    if( HermConj != 0){ HermConj = 0; }
    retlen = ( xlen / 2 ) + 1;
  }
  if( ret == R_NilValue || ret == R_MissingArg ){
    PROTECT(ret = Rf_allocVector(CPLXSXP, retlen));
    nprot++;
  } else {
    if( TYPEOF(ret) != CPLXSXP ){
      stop("ravetools `fftw_r2c`: `ret` should be complex");
    }
    if( Rf_xlength(ret) < retlen ){
      stop("ravetools `fftw_r2c`: `ret` length should be at least " + std::to_string(retlen));
    }
  }

  if( TYPEOF(data) != REALSXP ){
    PROTECT(data = Rf_coerceVector(data, REALSXP));
    nprot++;
  // } else if (MAYBE_REFERENCED(data)) {
  } else if(!inplace && fftwplanopt <= 0) {
    data = PROTECT(Rf_duplicate(data));
    nprot++;
  }

  cfft_r2c(&xlen, REAL(data), reinterpret_cast<fftw_complex*>(&COMPLEX(ret)[0]), &HermConj,
           &fftwplanopt);

  if(nprot > 0){
    UNPROTECT(nprot);
  }
  return ret;
}

// [[Rcpp::export]]
SEXP mvfftw_r2c(SEXP data,
               int fftwplanopt = 0,
               SEXP ret = R_NilValue,
               bool inplace = false)
{
  int nprot = 0;

  // check HermConj and ret
  int nrows = Rf_nrows(data);
  int ncols = Rf_ncols(data);
  int retrows = ( nrows / 2 ) + 1;
  if( ret == R_NilValue || ret == R_MissingArg ){
    PROTECT(ret = Rf_allocMatrix(CPLXSXP, retrows, ncols));
    nprot++;
  } else {
    if( TYPEOF(ret) != CPLXSXP ){
      stop("ravetools `fftw_r2c`: `ret` should be complex");
    }
    if( Rf_xlength(ret) != retrows * ncols ){
      stop("ravetools `fftw_r2c`: `ret` length should be " + std::to_string(retrows * ncols));
    }
  }



  if( TYPEOF(data) != REALSXP ){
    PROTECT(data = Rf_coerceVector(data, REALSXP));
    nprot++;
  // } else if (MAYBE_REFERENCED(data)) {
  } else if(!inplace && fftwplanopt <= 0) {
    // data need to be copied
    // however fftwplanopt > 0 will copy eventually, so only
    // copy when fftwplanopt <= 0
    data = PROTECT(Rf_duplicate(data));
    nprot++;
  }

  cmvfft_r2c(&nrows, &ncols, REAL(data),
             reinterpret_cast<fftw_complex*>(&COMPLEX(ret)[0]),
             &fftwplanopt);

  if(nprot > 0){
    UNPROTECT(nprot);
  }

  return(ret);
}

// [[Rcpp::export]]
SEXP fftw_c2c(SEXP data, int inverse = 0, SEXP ret = R_NilValue, bool inplace = false)
{
  int nprot = 0;
  int xlen = Rf_length(data);
  if(TYPEOF(data) != CPLXSXP){
    PROTECT(data = Rf_coerceVector(data, CPLXSXP));
    nprot++;
  // } else if (MAYBE_REFERENCED(data)) {
  } else if(!inplace) {
    data = PROTECT(Rf_duplicate(data));
    nprot++;
  }

  if(inverse){
    inverse = 1;
  }

  if(ret == R_NilValue){
    PROTECT(ret = Rf_allocVector(CPLXSXP, xlen));
    nprot++;
  } else {
    if(TYPEOF(ret) != CPLXSXP){
      stop("ravetools `fftw_c2c`: `ret` must be complex");
    }
    if(Rf_length(ret) != xlen) {
      stop("ravetools `fftw_c2c`: `ret` must have length of " + std::to_string(xlen));
    }
  }
  cfft_c2c(
    &xlen, reinterpret_cast<fftw_complex*>(&COMPLEX(data)[0]),
    reinterpret_cast<fftw_complex*>(&COMPLEX(ret)[0]), &inverse
  );

  if(nprot > 0){
    UNPROTECT(nprot);
  }

  return ret;
}

// [[Rcpp::export]]
SEXP fftw_c2r(SEXP data, int HermConj = 1, SEXP ret = R_NilValue, bool inplace = false){
  int nprot = 0;

  // check HermConj and ret
  int xlen = Rf_length(data);
  int retlen = 0;
  if( HermConj == 1 ){
    retlen = xlen;
  } else {
    if( HermConj != 0){ HermConj = 0; }
    retlen = (xlen - 1) * 2;
  }
  if( ret == R_NilValue || ret == R_MissingArg ){
    PROTECT(ret = Rf_allocVector(REALSXP, retlen));
    nprot++;
  } else {
    if( TYPEOF(ret) != REALSXP ){
      stop("ravetools `fftw_c2r`: `ret` should be double");
    }
    if( Rf_xlength(ret) < retlen ){
      stop("ravetools `fftw_c2r`: `ret` length should be at least " + std::to_string(retlen));
    }
    if( Rf_xlength(ret) > retlen ){
      retlen++;
    }
  }

  if( TYPEOF(data) != CPLXSXP ){
    PROTECT(data = Rf_coerceVector(data, CPLXSXP));
    nprot++;
  // } else if (MAYBE_REFERENCED(data)) {
  } else if(!inplace) {
    data = PROTECT(Rf_duplicate(data));
    nprot++;
  }

  cfft_c2r(&xlen, reinterpret_cast<fftw_complex*>(&COMPLEX(data)[0]),
           REAL(ret));

  if(nprot > 0){
    UNPROTECT(nprot);
  }
  return ret;
}

// [[Rcpp::export]]
SEXP conjugate(SEXP data) {
  if(TYPEOF(data) != CPLXSXP){
    stop("`conjugate`: data must be complex");
  }

  int xlen = Rf_length(data);
  int i = 0;
  for(Rcomplex* ptr = COMPLEX(data); i < xlen; i++, ptr++){
    ptr->i = -(ptr->i);
  }

  return R_NilValue;
}


/*** R
x <- rnorm(1000)
ret = double(1000)
a <- fftw_c2r(x, ret = ret)
b <- fftwtools::fftw_c2r(x, 1)
# max(Mod(b-Conj(a)))
range(b-a)
*/
