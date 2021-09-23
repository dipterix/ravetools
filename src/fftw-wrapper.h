#ifndef RAVEUTILS_FFTW_WRAPPER_H
#define RAVEUTILS_FFTW_WRAPPER_H

#include <Rcpp.h>

SEXP fftw_r2c(SEXP data, int HermConj, SEXP ret);

SEXP mvfft_r2c(SEXP data, int HermConj, int fftwplanopt, SEXP ret);

SEXP fftw_c2c(SEXP data, int inverse, SEXP ret);

SEXP conjugate(SEXP data);

#endif // RAVEUTILS_FFTW_WRAPPER_H
