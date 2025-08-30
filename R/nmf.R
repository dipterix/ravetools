#' @title A naive implementation of non-negative matrix factorization
#' @description
#' A pure-R vanilla implementation assuming inputs are non-negative matrices
#' without \code{NA}.
#' @param x a matrix, or can be converted into a matrix; all negative or missing
#' values will be treated as zero
#' @param k decomposition rank
#' @param tol stop criteria, a numeric of two; the first number is the
#' tolerance for root-mean-squared residuals, relative to the largest number in
#' \code{x}; the second number is the tolerance for weight differences; any
#' stopping criteria met will result in the stop of iteration
#' @param max_iters maximum iterations
#' @param verbose whether to report the progress; logical or a positive integer
#' (of step intervals)
#' @returns A list of weights (non-negative template matrix \code{W} and
#' non-negative \code{H}) and errors (root mean squared error of fitted,
#' matrix \code{W}, and \code{W} versus their previous iteration, respectively).
#'
#' @examples
#'
#'
#'
#' x <- stats::toeplitz(.9 ^ (0:31))
#'
#' nmf <- naive_nmf(x, k = 7, verbose = FALSE)
#'
#' fitted <- nmf$W %*% nmf$H
#'
#' oldpar <- par(mfrow = c(1, 2))
#' on.exit({ par(oldpar )})
#'
#' image(x, zlim = c(0, 1), main = "Input")
#' image(fitted, zlim = c(0, 1),
#'       main = sprintf("Fitted with rank=%d", nmf$k))
#'
#'
#'
#'
#' @export
naive_nmf <- function(x, k, tol = c(1e-4, 1e-8), max_iters = 1e4, verbose = TRUE) {
  max_iters <- as.integer(max_iters)
  stopifnot(isTRUE(max_iters >= 0))
  if(max_iters == 0) {
    max_iters <- 1e4
  }

  if(!is.matrix(x)) {
    x <- as.matrix(x)
  }
  dimnames(x) <- NULL
  x[is.na(x) | x < 0] <- 0
  tol1 <- max(x) * tol[[1]]
  tol2 <- tol[[2]]

  dm <- dim(x)
  nr <- dm[[1]]
  nc <- dm[[2]]

  results <- new.env(parent = emptyenv())
  results$k <- k
  results$iters <- 0L
  # nr x k
  results$W <- matrix(stats::runif(nr * k), nrow = nr, ncol = k)

  # k x nc
  results$H <- matrix(stats::runif(nc * k), nrow = k, ncol = nc)

  results$WH <- results$W %*% results$H

  results$error <- c(Inf, Inf, Inf)
  results$converged <- FALSE

  if(isTRUE(verbose)) {
    verbose_interval <- ceiling(max(max_iters / 100, 10))
    verbose <- TRUE
  } else if(!isFALSE(verbose)) {
    verbose_interval <- min(max(as.integer(verbose), 1), max_iters)
    verbose <- TRUE
  }

  eps <- .Machine$double.eps

  # errs <- dipsaus::fastqueue2()

  update <- function() {

    tmp <- crossprod(results$W, results$WH)
    diag(tmp) <- diag(tmp) + eps
    new_H <- results$H * ( crossprod(results$W, x) / tmp )

    new_H <- new_H * 0.1 + results$H * 0.9

    tmp <- tcrossprod( results$W %*% new_H, new_H )
    diag(tmp) <- diag(tmp) + eps
    new_W <- results$W * ( tcrossprod(x, new_H) / tmp )

    new_W <- new_W * 0.1 + results$W * 0.9

    new_WH <- new_W %*% new_H

    err_H <- sqrt(mean((new_H - results$H)^2, na.rm = FALSE))
    err_W <- sqrt(mean((new_W - results$W)^2, na.rm = FALSE))
    err_WH <- sqrt(mean((x - new_WH)^2, na.rm = FALSE))
    new_errors <- c(err_WH, err_W, err_H)

    results$H <- new_H
    results$W <- new_W
    results$WH <- new_WH

    # results$error_diff <- new_errors - results$error
    results$error <- new_errors
    # errs$add(new_errors)

    results$iters <- results$iters + 1L
    return()
  }

  lapply(seq_len(max_iters), function(iter) {
    if( results$converged ) { return() }
    update()

    new_errors <- results$error

    if((new_errors[[2]] < tol2 && new_errors[[3]] < tol2) || new_errors[[1]] < tol1) {
      results$converged <- TRUE
    }

    if( verbose && (iter %% verbose_interval == 0 || iter == max_iters || results$converged) ) {
      cat(sprintf("Iteration %04d: [fitted err=%.2e] [W err=%.2e] [H err=%.2e]        \r",
                  iter, new_errors[[1]], new_errors[[2]], new_errors[[3]]))
    }
    return()

  })
  if(verbose) {
    cat("\n")
  }

  # results$error_trace <- do.call("rbind", errs$as_list())
  as.list(results)

}
