#include "dijkstras-path.h"
#include "TinyParallel.h"
using namespace Rcpp;

// assuming they have the same length
double pointDistance(const double* a, const double *b, const size_t len) {
  double re = 0.0;
  // Rcout << len << "!\n";
  for(size_t i = 0; i < len; i++) {
    // Rcout << *(a + i) << " " << *(b + i) << "\n";
    re += std::pow((*(a + i) - *(b + i)), 2.0);
  }
  return std::sqrt(re);
}

// [[Rcpp::export]]
SEXP dijkstras_path(
    const SEXP& position,
    const SEXP& index,
    const size_t& nPoints,
    const size_t& nIndices,
    const size_t& startIndex,
    const double& maxDistance,
    const double& maxEdgeLen
) {
  SEXP re = R_NilValue;

  bool maxEdgeLenUnset = maxEdgeLen <= 0.0;

  // validate data
  const R_xlen_t posLen = Rf_xlength(position);
  size_t pointSize = ((size_t) posLen) / nPoints;
  if ( pointSize <= 0 || pointSize * nPoints - posLen != 0 ) {
    re = PROTECT(make_error("C++ `dijkstras_path`: `position` length is not a multiple of `nPoints`"));
    UNPROTECT(1);
    return re;
  } else {
    Rcout << "Detected point size: " << pointSize << "\n";
  }
  const R_xlen_t idxLen = Rf_xlength(index);
  const R_xlen_t faceSize = idxLen / (R_xlen_t) nIndices;
  if ( faceSize < 1 || faceSize * nIndices - idxLen != 0 ) {
    re = PROTECT(make_error("C++ `dijkstras_path`: `index` length is not a multiple of `nIndices`, or `faceSize` is too small (needs to be at least 1)"));
    UNPROTECT(1);
    return re;
  }

  SEXP position_;
  SEXP index_;
  if( TYPEOF(position) - REALSXP != 0 ) {
    position_ = PROTECT(Rf_coerceVector(position, REALSXP));
  } else {
    position_ = PROTECT(position);
  }
  if( TYPEOF(index) - INTSXP != 0 ) {
    index_ = PROTECT(Rf_coerceVector(index, INTSXP));
  } else {
    index_ = PROTECT(index);
  }
  int* ptrIndex = INTEGER(index_);
  double* ptrPosition = REAL(position_);

  // position is pointSize x nPoints matrix
  // index is faceSize x idxLen matrix

  // prepare table and data
  double* ptrCurrent = ptrPosition + (startIndex * pointSize);
  R_xlen_t idxCurrent = (R_xlen_t) startIndex;

  // used as temporary iterator index
  R_xlen_t ii, idxTmp;
  R_xlen_t iiMask = faceSize - 1;
  double edgeDistance, minDistance, currentDistance;

  // table columns are: node (target, automatically from 0 - posLen-1),
  // previous node, distance (from target to previous node), visited
  // (whether the distance is complete)
  SEXP prevNode = PROTECT( Rf_allocVector( INTSXP, nPoints ) );
  SEXP distance = PROTECT( Rf_allocVector( REALSXP, nPoints ) );
  SEXP visited = PROTECT( Rf_allocVector( RAWSXP, nPoints ) );

  re = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(re, 0, prevNode);
  SET_VECTOR_ELT(re, 1, distance);

  // initialize prevNode with NA
  int* ptrPrevNode = INTEGER(prevNode);
  double* ptrDistance0 = REAL(distance);
  double* ptrDistance = REAL(distance);
  Rbyte* ptrVisited0 = RAW(visited);
  Rbyte* ptrVisited = RAW(visited);

  for( ii = 0; ii < nPoints; ii++ ) {
    *ptrPrevNode++ = NA_INTEGER;
    *ptrDistance++ = -1.0;
    *ptrVisited++ = 0;
  }

  // reset and reuse these pointers
  ptrPrevNode = INTEGER(prevNode);
  ptrDistance = ptrDistance0;
  ptrVisited = ptrVisited0;

  // initialize the first node, which is the starting node
  *(ptrDistance0 + idxCurrent) = 0.0;
  *(ptrVisited0 + idxCurrent) = 1;


  // loop
  for(size_t loopIdx = 0; loopIdx < nPoints; loopIdx++) {

    if( (loopIdx * 1023) == 0 ) {
      Rcpp::checkUserInterrupt();
    }

    // Find node visited (not finished) with the shortest distance
    minDistance = -1.0;
    // Rcpp::print(visited);
    for( ii = 0, ptrVisited = ptrVisited0, ptrDistance = ptrDistance0; ii < nPoints; ii++, ptrVisited++, ptrDistance++ ) {
      if( *ptrVisited == 1 && (minDistance > *ptrDistance || minDistance < 0.0) ) {
        idxCurrent = ii;
        minDistance = *ptrDistance;
        // Rcout << "Current node is: " << idxCurrent << " (dist=" << minDistance << ")"<< "\n";
      }
    }

    if( minDistance < 0.0 || (maxDistance > 0.0 && minDistance >= maxDistance) ) {
      break;
    }

    // reset face index pointer
    ptrIndex = INTEGER(index_);

    // minDistance = -1.0;

    // current distance to starting point
    currentDistance = *( ptrDistance0 + idxCurrent );
    ptrCurrent = ptrPosition + (idxCurrent * pointSize);
    Rcout << loopIdx << ": Current node is: " << idxCurrent << " (dist=" << currentDistance << ")       \n";

    for(ii = 0; ii < idxLen; ii++, ptrIndex++) {
      // find adjacent node with NA prevNode
      if( *ptrIndex == idxCurrent ) {
        if( (ii & iiMask) > 0 ) {
          idxTmp = *(ptrIndex - 1);
          ptrVisited = ptrVisited0 + idxTmp;
          if( *ptrVisited < 2 ) {
            edgeDistance = pointDistance(ptrCurrent, ptrPosition + (idxTmp * pointSize), pointSize);
            ptrDistance = ptrDistance0 + idxTmp;

            if( maxEdgeLenUnset || edgeDistance < maxEdgeLen ) {
              if( *ptrDistance < 0.0 || *ptrDistance > (currentDistance + edgeDistance) ) {
                // Rcout << "Updating node: " << idxTmp << " " <<
                //   *( ptrPrevNode + idxTmp ) << " (dist=" << *ptrDistance << ") -> " <<
                //     idxCurrent << " (dist=" << currentDistance << " + "<< edgeDistance << ")\n";
                *ptrDistance = currentDistance + edgeDistance;
                *( ptrPrevNode + idxTmp ) = idxCurrent;
              }
              // if( minDistance > edgeDistance || minDistance < 0.0 ) {
              //   minDistance = edgeDistance;
              // }
              *ptrVisited = 1;
            }
          }
        }
        if( (ii & iiMask) < iiMask ) {
          idxTmp = *(ptrIndex + 1);
          ptrVisited = ptrVisited0 + idxTmp;
          if( *ptrVisited < 2 ) {
            edgeDistance = pointDistance(ptrCurrent, ptrPosition + (idxTmp * pointSize), pointSize);
            ptrDistance = ptrDistance0 + idxTmp;
            if( maxEdgeLenUnset || edgeDistance < maxEdgeLen ) {
              if( *ptrDistance < 0.0 || *ptrDistance > (currentDistance + edgeDistance) ) {
                // Rcout << "Updating node: " << idxTmp << " " <<
                //   *( ptrPrevNode + idxTmp ) << " (dist=" << *ptrDistance << ") -> " <<
                //     idxCurrent << " (dist=" << currentDistance << " + "<< edgeDistance << ")\n";
                *ptrDistance = currentDistance + edgeDistance;
                *( ptrPrevNode + idxTmp ) = idxCurrent;
              }
              // if( minDistance > edgeDistance || minDistance < 0.0 ) {
              //   minDistance = edgeDistance;
              // }
              *ptrVisited = 1;
              // Rcout << "Current node: " << idxCurrent << " Found node " << idxTmp << " (dist=" << *ptrDistance << ")\n";
            }
          }
        }
      }
    }

    *(ptrVisited0 + idxCurrent) = 2;
  }

  UNPROTECT(6); // re, visited, distance, prevNode, index_, position_

  return re;
}

/*** R
devtools::load_all()
position <- matrix(nrow = 2, byrow = FALSE, data = c(
  0.0, 0.0,
  1, 0,
  0, 2,
  1, 1
))
index <- matrix(nrow = 2, byrow = FALSE, data = c(
  0, 2,
  0, 3,
  1, 2,
  1, 3
))
dijkstras_path(
  position, index = index[c(2,1), ], nPoints = 4L, nIndices = 4L, startIndex = 0L
)
mesh <- freesurferformats::read.fs.surface("~/rave_data/others/three_brain/N27/surf/lh.pial")
re <- dijkstras_path(
  position = t(mesh$vertices),
  index = t(mesh$faces) - min(mesh$faces),
  nPoints = nrow(mesh$vertices),
  nIndices = nrow(mesh$faces),
  startIndex = 0L,
  maxDistance = 521
)

*/
