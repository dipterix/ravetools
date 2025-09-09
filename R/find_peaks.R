#' Find peaks of a signal
#' @param x a numeric vector without missing values
#' @param min_val find peaks that are greater than this value
#' @param min_distance merge peaks that are less than \code{min_distance}
#' time-points away
#' @param min_width search radius (time-points) on whether the peak is
#' "local"; this is for seasonal oscillations.
#' @returns A list of peak index (1-based) and the corresponding value.
#' @examples
#'
#'
#' # Basic example
#' x <- sin(seq(0, 10, 0.01)) + rnorm(1001) * 0.1
#'
#' peaks <- find_peaks(x)
#'
#' plot(x, type = 'l')
#' abline(v = peaks$index, col = 'red')
#'
#' # merge peaks that are close
#' peaks <- find_peaks(x, min_distance = 400)
#'
#' plot(x, type = 'l')
#' abline(v = peaks$index, col = 'red')
#'
#' # with or without min_width
#' x <- c(0, 1, 0.5, 0.9, 0.2, 0.8, 0.2, 0.75, 0)
#'
#' # without min_width
#' peaks <- find_peaks(x, min_width = 0)
#' plot(x, type = 'l')
#' abline(v = peaks$index, col = 'red')
#'
#'
#' # with min_width=2: t=4 is greater than t=6
#' peaks <- find_peaks(x, min_width = 2)
#' plot(x, type = 'l')
#' abline(v = peaks$index, col = 'red')
#'
#'
#' @export
find_peaks <- function(x, min_val = NA, min_distance = 4, min_width = 2){
  df1 <- c(1, diff(x))
  df2 <- c(1, 1, diff(x, differences = 2))

  # necessary condition to be peaks (1st & 2nd derivatives)
  idx <- (df1 * c(df1[-1], 0)) <= 0 & c(df2[-1], 0) < 0

  # Get peaks that are beyond min_val
  if( is.na(min_val) ){
    # automatically decide min_val from standard error
    min_val <- sd(x)
    if(all(x >= 0)){
      min_val <- min(x) + min_val
    } else {
      min_val <- median(x) + min_val
    }
  }
  idx <- which(idx & (x > min_val))

  ord <- order(x[idx], decreasing = TRUE)
  idx_desc <- idx[ord]

  # merge peaks that are less than min_distance away
  for(ii in seq_along(idx_desc)){
    elem <- idx_desc[[ii]]
    idx_desc[idx_desc > 0 & abs(idx_desc - elem) < min_distance] <- elem
  }
  idx_desc <- unique(idx_desc)

  # check +- min_width to see if the it is still a peak
  p_idx <- idx_desc + min_width
  p_idx[p_idx > length(x)] <- length(x)
  m_idx <- idx_desc - min_width
  m_idx[m_idx < 1] <- 1
  sign <- (x[p_idx] - x[idx_desc]) * (x[idx_desc] - x[m_idx])
  idx_desc <- idx_desc[sign <= 0]
  list(
    index = idx_desc,
    values = x[idx_desc]
  )
}
