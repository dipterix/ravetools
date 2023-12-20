#ifndef RAVETOOLS_DIJKSTRA_H
#define RAVETOOLS_DIJKSTRA_H

#include "utils.h"
SEXP dijkstras_path(
    const SEXP& position,
    const SEXP& index,
    const SEXP& indexOrder,
    const size_t& nPoints,
    const size_t& nIndices,
    const size_t& startIndex,
    const double& maxDistance = 0.0,
    const double& maxEdgeLen = 0.0,
    const bool& verbose = true
);

#endif // RAVETOOLS_DIJKSTRA_H
