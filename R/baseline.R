#' @title Calculate Contrasts of Arrays in Different Methods
#' @description Provides five methods to baseline an array and calculate
#' contrast.
#' @param x array (tensor) to calculate contrast
#' @param along_dim integer range from 1 to the maximum dimension of \code{x}.
#' baseline along this dimension, this is usually the time dimension.
#' @param baseline_indexpoints integer vector, which index points are counted
#' into baseline window? Each index ranges from 1 to \code{dim(x)[[along_dim]]}.
#' See Details.
#' @param unit_dims integer vector, baseline unit: see Details.
#' @param method character, baseline method options are:
#' \code{"percentage"}, \code{"sqrt_percentage"}, \code{"decibel"},
#' \code{"zscore"}, and \code{"sqrt_zscore"}
#' @param baseline_subarray sub-arrays that should be used to calculate
#' baseline; default is \code{NULL} (automatically determined by
#' \code{baseline_indexpoints}).
#' @param ... passed to other methods
#'
#' @returns Contrast array with the same dimension as \code{x}.
#' @details
#' Consider a scenario where we want to baseline a bunch of signals recorded
#' from different locations. For each location, we record \code{n} sessions.
#' For each session, the signal is further decomposed into frequency-time
#' domain. In this case, we have the input \code{x} in the following form:
#' \deqn{session x frequency x time x location}
#' Now we want to calibrate signals for each session, frequency and location
#' using the first 100 time points as baseline points, then the code will be
#' \deqn{baseline_array(x, along_dim=3, baseline_window=1:100, unit_dims=c(1,2,4))}
#' \code{along_dim=3} is dimension of time, in this case, it's the
#' third dimension of \code{x}. \code{baseline_indexpoints=1:100}, meaning
#' the first 100 time points are used to calculate baseline.
#' \code{unit_dims} defines the unit signal. Its value \code{c(1,2,4)}
#' means the unit signal is per session (first dimension), per frequency
#' (second) and per location (fourth).
#'
#' In some other cases, we might want to calculate baseline across frequencies
#' then the unit signal is \eqn{frequency x time}, i.e. signals that share the
#' same session and location also share the same baseline. In this case,
#' we assign \code{unit_dims=c(1,4)}.
#'
#' There are five baseline methods. They fit for different types of data.
#' Denote \eqn{z} is an unit signal, \eqn{z_0} is its baseline slice. Then
#' these baseline methods are:
#'
#' \describe{
#' \item{\code{"percentage"}}{
#' \deqn{
#'   \frac{z - \bar{z_{0}}}{\bar{z_{0}}} \times 100\%
#' }{
#'   (z / mean(z_0) - 1) x 100\%
#' }
#' }
#' \item{\code{"sqrt_percentage"}}{
#' \deqn{
#'   \frac{\sqrt{z} - \bar{\sqrt{z_{0}}}}{\bar{\sqrt{z_{0}}}} \times 100\%
#' }{
#'   (sqrt(z) / mean(sqrt(z_0)) - 1) x 100\%
#' }
#' }
#' \item{\code{"decibel"}}{
#' \deqn{
#'   10 \times ( \log_{10}(z) - \bar{\log_{10}(z_{0})} )
#' }{
#'   10 * ( log10 (z) - mean( log10(z_0) ) )
#' }
#' }
#' \item{\code{"zscore"}}{
#' \deqn{
#'   \frac{z-\bar{z_{0}}}{sd(z_{0})}
#' }{
#'   (z - mean( z_0 )) / sd( z_0 )
#' }
#' }
#' \item{\code{"sqrt_zscore"}}{
#' \deqn{
#'   \frac{\sqrt{z}-\bar{\sqrt{z_{0}}}}{sd(\sqrt{z_{0}})}
#' }{
#'   (sqrt(z) - mean( sqrt(z_0) )) / sd( sqrt(z_0) )
#' }
#' }
#'
#'
#' }
#'
#' @examples
#'
#' # Set ncores = 2 to comply to CRAN policy. Please don't run this line
#' ravetools_threads(n_threads = 2L)
#'
#'
#' library(ravetools)
#' set.seed(1)
#'
#' # Generate sample data
#' dims = c(10,20,30,2)
#' x = array(rnorm(prod(dims))^2, dims)
#'
#' # Set baseline window to be arbitrary 10 timepoints
#' baseline_window = sample(30, 10)
#'
#' # ----- baseline percentage change ------
#'
#' # Using base functions
#' re1 <- aperm(apply(x, c(1,2,4), function(y){
#'   m <- mean(y[baseline_window])
#'   (y/m - 1) * 100
#' }), c(2,3,1,4))
#'
#' # Using ravetools
#' re2 <- baseline_array(x, 3, c(1,2,4),
#'                       baseline_indexpoints = baseline_window,
#'                       method = 'percentage')
#'
#' # Check different, should be very tiny (double precisions)
#' range(re2 - re1)
#'
#' \donttest{
#' # Check speed for large dataset, might take a while to profile
#'
#' ravetools_threads(n_threads = -1)
#'
#' dims <- c(200,20,300,2)
#' x <- array(rnorm(prod(dims))^2, dims)
#' # Set baseline window to be arbitrary 10 timepoints
#' baseline_window <- seq_len(100)
#' f1 <- function(){
#'   aperm(apply(x, c(1,2,4), function(y){
#'     m <- mean(y[baseline_window])
#'     (y/m - 1) * 100
#'   }), c(2,3,1,4))
#' }
#' f2 <- function(){
#'   # equivalent as bl = x[,,baseline_window, ]
#'   #
#'   baseline_array(x, along_dim = 3,
#'                  baseline_indexpoints = baseline_window,
#'                  unit_dims = c(1,2,4), method = 'percentage')
#' }
#' range(f1() - f2())
#' microbenchmark::microbenchmark(f1(), f2(), times = 10L)
#'
#' }
#'
#'
#'
#' @export
baseline_array <- function(x, along_dim, unit_dims = seq_along(dim(x))[-along_dim], ...) {

  dims <- dim(x)
  along_dim <- as.integer(along_dim)
  unit_dims <- as.integer(unit_dims)

  stopifnot2(along_dim >=1 && along_dim <= length(dims), msg = paste0(
    sQuote('along_dim'), ' is invalid, it must be an integer from 1 to ',
    length(dims)
  ))

  # unit_dims is baseline unit,
  stopifnot2(!any(is.na(unit_dims)), msg = paste0(sQuote('unit_dims'), ' contains NAs'))
  stopifnot2(all(unit_dims %in% seq_along(dims)), msg = paste0(sQuote('unit_dims'), ' has invalid dimensions'))
  stopifnot2(!along_dim %in% unit_dims, msg = paste0(
    sQuote('along_dim'), ' cannot be inside of ', sQuote('unit_dims')))
  unit_dims <- sort(unit_dims)

  UseMethod("baseline_array")

}

#' @rdname baseline_array
#' @export
baseline_array.array <- function(
    x, along_dim, unit_dims = seq_along(dim(x))[-along_dim],
    method = c("percentage", "sqrt_percentage", "decibel", "zscore",
               "sqrt_zscore", "subtract_mean"),
    baseline_indexpoints = NULL, baseline_subarray = NULL, ...) {

  along_dim <- as.integer(along_dim)
  unit_dims <- as.integer(unit_dims)
  method <- match.arg(method)
  method_int <- which(c(
    "percentage", "sqrt_percentage", "decibel", "zscore", "sqrt_zscore",
    "subtract_mean") == method)

  dims <- dim(x)
  ntimepoints <- dims[along_dim]

  stopifnot2(!(is.null(baseline_indexpoints) && is.null(baseline_subarray)),
             msg = "Either `baseline_indexpoints` or `baseline_subarray` must be specified")
  if(!is.null(baseline_indexpoints)) {
    baseline_indexpoints <- as.integer(baseline_indexpoints)
    # calculate baseline window
    baseline_indexpoints <- baseline_indexpoints[!is.na(baseline_indexpoints) & baseline_indexpoints > 0 & baseline_indexpoints <= ntimepoints]
    baseline_indexpoints <- sort(baseline_indexpoints)
    stopifnot2(length(baseline_indexpoints) > 0, msg = paste0('Baseline window is invalid: cannot find any valid time points. \nPlease makesure ', sQuote('baseline_indexpoints'), ' contains at least one integer from 1 to ', ntimepoints))
    call_args <- lapply(seq_along(dims), function(ii){
      if(ii == along_dim){
        return(quote(baseline_indexpoints))
      }
      .missing_arg[[1]]
    })
    call_args <- c(list(quote(x)), call_args, list(drop = FALSE))
    bl <- do.call('[', call_args)
    bldims <- dim(bl)
  } else {
    stopifnot2(
      is.numeric(baseline_subarray) && is.array(baseline_subarray) &&
        length(dim(baseline_subarray)) == length(dims),
      msg = "`baseline_subarray` must be an array with same margin size as `x`"
    )
    bldims2 <- dim(baseline_subarray)
    bldims <- dims
    bldims[[along_dim]] <- bldims2[[along_dim]]
    bldims <- as.integer(bldims)
    stopifnot2(all(bldims2 == bldims),
               msg = "`baseline_subarray` dimension is invalid")
    bl <- baseline_subarray
  }

  rest <- seq_along(dims)[-unit_dims]
  #
  # if(!all(c(1, 3, 4) %in% unit_dims)) {
  #   bl <- ravetools::collapse(bl, keep = unit_dims)
  # }
  #
  #
  # dim(bl2) <- c(20,1,1,1)
  # bldims2 <- dim(bl2)
  # re <- .Call(`_ravetools_baselineArray`, x, bl2, dims, bldims2, along_dim - 1L, unit_dims - 1L, rest - 1L, method_int - 1L)

  re <- .Call(`_ravetools_baselineArray`, x, bl, dims, bldims, along_dim - 1L, unit_dims - 1L, rest - 1L, method_int - 1L)
  if(inherits(re, "ravetools_error")){
    stop(re)
  }
  re
}
