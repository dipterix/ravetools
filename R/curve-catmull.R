# --------------------------------------------------------------------------- #
#  Catmull-Rom 3D spline – mirrors THREE.CatmullRomCurve3
#  Supports: centripetal (default), chordal, uniform (catmullrom + tension)
# --------------------------------------------------------------------------- #

# ---- cubic Hermite polynomial helpers (internal) --------------------------

# Coefficients of  f(t) = c0 + c1*t + c2*t^2 + c3*t^3
# with f(0)=x0, f(1)=x1, f'(0)=v0, f'(1)=v1.
crm_cubic_init <- function(x0, x1, v0, v1) {
  c(x0,
    v0,
    -3 * x0 + 3 * x1 - 2 * v0 - v1,
     2 * x0 - 2 * x1 +     v0 + v1)
}

crm_cubic_eval <- function(coef, t) {
  coef[[1L]] + t * (coef[[2L]] + t * (coef[[3L]] + t * coef[[4L]]))
}

# Non-uniform Barry-Goldman reparametrization
crm_nonuniform <- function(x0, x1, x2, x3, dt0, dt1, dt2) {
  v1 <- (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1
  v2 <- (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2
  crm_cubic_init(x1, x2, v1 * dt1, v2 * dt1)
}

# Classic uniform Catmull-Rom (tension-scaled)
crm_uniform <- function(x0, x1, x2, x3, tension) {
  crm_cubic_init(x1, x2, tension * (x2 - x0), tension * (x3 - x1))
}

# --------------------------------------------------------------------------- #

#' @title \verb{Catmull-Rom} 3D Spline Curve
#' @description
#' Creates a smooth \verb{Catmull-Rom} spline curve through a set of 3D key points.
#'
#' @param points numeric matrix with at least 2 rows and exactly 3 columns
#'   (\code{x}, \code{y}, \code{z}), giving the key (control) points through
#'   which the curve passes.
#' @param curve_type character; One of
#'   \code{"centripetal"} (default), \code{"chordal"}, or \code{"uniform"}.
#'   \code{"centripetal"} uses \eqn{\alpha = 0.5} (square-root of chord
#'   length), \code{"chordal"} uses \eqn{\alpha = 1} (chord length), and
#'   \code{"uniform"} is the classic formulation controlled by
#'   \code{tension}.
#' @param tension numeric scalar in \eqn{[0, 1]}; tangent scaling factor used
#'   only when \code{curve_type = "uniform"}. At \code{0.5} (default) the
#'   curve matches the standard \verb{Catmull-Rom} formulation.
#' @param closed logical; if \code{TRUE} the curve closes on itself by
#'   connecting the last point back to the first. Default is \code{FALSE}.
#'
#' @returns An object of class \code{"ravetools_curve"} (a list) with the
#'   following elements:
#' \describe{
#'   \item{\code{points}}{The input key-point matrix (\eqn{n \times 3}).}
#'   \item{\code{curve_type}}{Character, the parameterization type.}
#'   \item{\code{tension}}{Numeric, the tension value (relevant for
#'     \code{"uniform"} only).}
#'   \item{\code{closed}}{Logical, whether the curve is closed.}
#'   \item{\code{get_point}}{A \code{function(t)} that accepts a scalar
#'     \code{t} in \eqn{[0, 1]} and returns a named numeric vector
#'     on the curve.}
#'   \item{\code{get_points}}{A \code{function(n)} that returns an
#'     \eqn{n \times 3} matrix of \code{n} evenly spaced points along the
#'     curve, with column names \code{"x"}, \code{"y"}, \code{"z"}.}
#'   \item{\code{t_keypoints}}{Numeric vector of length \eqn{n} with the
#'     \code{t} parameter value where each key point lies on the curve.
#'     First element is always \code{0}, last is always \code{1}.}
#'   \item{\code{segment_lengths}}{Numeric vector of length \eqn{n-1} (open
#'     curve) or \eqn{n} (closed curve) containing the arc length of each
#'     spline segment, estimated by numerical integration.}
#' }
#'
#' @seealso \code{\link{print.ravetools_curve}},
#'   \code{\link{plot.ravetools_curve}}
#'
#' @examples
#'
#' pts <- matrix(c(
#'   -33.0534, -10.6213, -21.8328,
#'   -34.7526, -25.5089, -14.5390,
#'   -41.2002, -10.4606, -22.0032,
#'   -46.4717, -10.3567, -22.1134,
#'   -51.7431, -10.2528, -22.2237,
#'   -57.0146, -10.1488, -22.3339,
#'   -62.2860, -10.0449, -22.4442,
#'   -67.5575,  -9.9410, -22.5544
#' ), ncol = 3, byrow = TRUE)
#'
#' curve <- catmull_rom_3d(pts)
#' print(curve)
#'
#' # Sample 100 evenly spaced points along the curve
#' smooth <- curve$get_points(100)
#' head(smooth)
#'
#' # Evaluate the curve at t = 0.5 (midpoint)
#' curve$get_point(0.5)
#'
#' plot(curve, use_rgl = FALSE)
#'
#' @export
catmull_rom_3d <- function(
    points,
    curve_type = c("centripetal", "chordal", "uniform"),
    tension    = 0.5,
    closed     = FALSE
) {
  curve_type <- match.arg(curve_type)

  if (!is.matrix(points)) {
    points <- as.matrix(points)
  }
  stopifnot2(is.numeric(points),
             msg = "`points` must be a numeric matrix.")
  stopifnot2(ncol(points) == 3L,
             msg = "`points` must have exactly 3 columns (x, y, z).")
  n_pts <- nrow(points)
  stopifnot2(n_pts >= 2L,
             msg = "`points` must have at least 2 rows.")

  tension <- as.numeric(tension)[[1L]]
  stopifnot2(is.finite(tension) && tension >= 0 && tension <= 1,
             msg = "`tension` must be a single finite number in [0, 1].")

  closed <- isTRUE(closed)

  # Number of key points captured at construction time
  l <- n_pts

  # t parameter values at each key point
  # Open:   t_i = i / (l-1)  => t[1]=0, t[l]=1
  # Closed: t_i = i / l      => t[1]=0, t[l]=(l-1)/l
  t_keypoints <- seq.int(0L, l - 1L) / (if (closed) l else (l - 1L))

  # Helper: 0-based index i → 1-based row, wrapping with modulo
  get_ctrl <- function(i) {
    points[(i %% l) + 1L, ]
  }

  # ---- get_point -----------------------------------------------------------
  get_point <- function(t) {
    t <- as.numeric(t[[1L]])

    # Map t ∈ [0,1] to floating segment index + intra-segment weight
    p_raw <- (l - (if (closed) 0L else 1L)) * t
    seg   <- floor(p_raw)
    w     <- p_raw - seg

    # Open curve: clamp t == 1 exactly onto the last segment
    if (!closed && w == 0 && seg == l - 1L) {
      seg <- l - 2L
      w   <- 1
    }

    # Four bounding control points (ghost extrapolation at open endpoints)
    # Note: R's %% is always non-negative for positive modulus, so negative
    # indices wrap correctly without extra adjustment (unlike JS).
    if (closed || seg > 0L) {
      p0 <- get_ctrl(seg - 1L)
    } else {
      # Reflect: p0 = 2*p1 - p2
      p0 <- 2 * points[1L, ] - points[2L, ]
    }

    p1 <- get_ctrl(seg)
    p2 <- get_ctrl(seg + 1L)

    if (closed || (seg + 2L) < l) {
      p3 <- get_ctrl(seg + 2L)
    } else {
      # Reflect: p3 = 2*p_{n} - p_{n-1}
      p3 <- 2 * points[l, ] - points[l - 1L, ]
    }

    # Build cubic Hermite polynomials per axis
    if (curve_type == "centripetal" || curve_type == "chordal") {
      # pow applied to squared-distance:
      #   centripetal: (d^2)^0.25 = d^0.5  (alpha = 0.5)
      #   chordal    : (d^2)^0.5  = d      (alpha = 1)
      pow_exp <- if (curve_type == "chordal") 0.5 else 0.25

      dt0 <- sum((p1 - p0)^2)^pow_exp
      dt1 <- sum((p2 - p1)^2)^pow_exp
      dt2 <- sum((p3 - p2)^2)^pow_exp

      # Safety check for degenerate (duplicate / very close) points
      if (dt1 < 1e-4) dt1 <- 1.0
      if (dt0 < 1e-4) dt0 <- dt1
      if (dt2 < 1e-4) dt2 <- dt1

      cx <- crm_nonuniform(p0[[1L]], p1[[1L]], p2[[1L]], p3[[1L]], dt0, dt1, dt2)
      cy <- crm_nonuniform(p0[[2L]], p1[[2L]], p2[[2L]], p3[[2L]], dt0, dt1, dt2)
      cz <- crm_nonuniform(p0[[3L]], p1[[3L]], p2[[3L]], p3[[3L]], dt0, dt1, dt2)
    } else {
      # uniform (classic Catmull-Rom with tension)
      cx <- crm_uniform(p0[[1L]], p1[[1L]], p2[[1L]], p3[[1L]], tension)
      cy <- crm_uniform(p0[[2L]], p1[[2L]], p2[[2L]], p3[[2L]], tension)
      cz <- crm_uniform(p0[[3L]], p1[[3L]], p2[[3L]], p3[[3L]], tension)
    }

    c(x = crm_cubic_eval(cx, w),
      y = crm_cubic_eval(cy, w),
      z = crm_cubic_eval(cz, w))
  }

  # ---- get_points ----------------------------------------------------------
  get_points <- function(n = 50L) {
    n <- as.integer(n[[1L]])
    if (is.na(n) || n < 2L) {
      stop("`n` must be an integer >= 2.")
    }
    ts  <- seq.int(0L, n - 1L) / (n - 1L)
    mat <- matrix(0, nrow = n, ncol = 3L)
    for (i in seq_len(n)) {
      mat[i, ] <- get_point(ts[[i]])
    }
    colnames(mat) <- c("x", "y", "z")
    mat
  }

  # ---- segment arc lengths (numerical integration, 50 sub-samples each) ---
  n_segs    <- if (closed) l else l - 1L
  t_seg_end <- c(t_keypoints[-1L], if (closed) 1.0)
  sub_n     <- 50L
  segment_lengths <- numeric(n_segs)
  for (s in seq_len(n_segs)) {
    ts_sub  <- seq(t_keypoints[[s]], t_seg_end[[s]], length.out = sub_n + 1L)
    pts_sub <- matrix(0, nrow = sub_n + 1L, ncol = 3L)
    for (j in seq_len(sub_n + 1L)) {
      pts_sub[j, ] <- get_point(ts_sub[[j]])
    }
    diffs <- diff(pts_sub)
    segment_lengths[[s]] <- sum(sqrt(rowSums(diffs * diffs)))
  }

  structure(
    list(
      points          = points,
      curve_type      = curve_type,
      tension         = tension,
      closed          = closed,
      t_keypoints     = t_keypoints,
      segment_lengths = segment_lengths,
      get_point       = get_point,
      get_points      = get_points
    ),
    class = "ravetools_curve"
  )
}

# --------------------------------------------------------------------------- #
#  S3 methods
# --------------------------------------------------------------------------- #

#' @title Print method for \code{ravetools_curve}
#' @description Prints a concise summary of a \code{ravetools_curve} object
#' returned by \code{\link{catmull_rom_3d}}.
#' @param x an object of class \code{ravetools_curve}.
#' @param ... additional arguments passed to \code{\link{print.data.frame}}.
#' @return Invisibly returns \code{x}.
#' @export
print.ravetools_curve <- function(x, ...) {
  type_str <- x$curve_type
  if (x$curve_type == "uniform") {
    type_str <- sprintf("uniform (tension=%.2f)", x$tension)
  }
  cat(sprintf("<ravetools_curve: %d key points, type=%s%s>\n",
              nrow(x$points),
              type_str,
              if (x$closed) ", closed" else ""))
  cat("Key points (t / x / y / z / seg_length):\n")
  seg_len <- c(x$segment_lengths, NA_real_)
  print(data.frame(
    t          = x$t_keypoints,
    x          = x$points[, 1L],
    y          = x$points[, 2L],
    z          = x$points[, 3L],
    seg_length = seg_len
  ), ...)
  cat(sprintf("Total arc length: %.4g\n", sum(x$segment_lengths)))
  invisible(x)
}

#' @title Plot method for \code{ravetools_curve}
#' @description
#' Plots a \code{ravetools_curve} object created by
#' \code{\link{catmull_rom_3d}}.  When the \code{rgl} package is available
#' and \code{use_rgl = TRUE} (default), an interactive 3D scene is opened.
#' Otherwise three 2-D projection panels (\emph{x-y}, \emph{x-z},
#' \emph{y-z}) are drawn using base \R{} graphics.
#'
#' @param x an object of class \code{ravetools_curve}.
#' @param n integer; number of sample points used to draw the smooth curve.
#'   Default is \code{200L}.
#' @param col color for the spline curve line. Default \code{"steelblue"}.
#' @param pch,cex plotting character and scaling for the key control points
#'   (base-R fallback only). Default \code{pch = 19}, \code{cex = 1}.
#' @param use_rgl logical; if \code{TRUE} (default) and the \code{rgl}
#'   package is installed, an interactive 3D window is used. Set to
#'   \code{FALSE} to force the base-R 2D projection panels.
#' @param ... additional graphical parameters forwarded to the underlying
#'   plot calls.
#' @return Invisibly returns \code{x}.
#' @export
plot.ravetools_curve <- function(x, n = 200L, col = "steelblue",
                                 pch = 19L, cex = 1,
                                 use_rgl = TRUE, ...) {
  smooth <- x$get_points(as.integer(n))
  kp     <- x$points

  # ---- rgl 3D path --------------------------------------------------------
  if (isTRUE(use_rgl) && check_rgl(strict = FALSE)) {
    return(rgl_view({
      rgl_call("lines3d",
               smooth[, 1L], smooth[, 2L], smooth[, 3L],
               col = col, lwd = 2, ...)
      rgl_call("points3d",
               kp[, 1L], kp[, 2L], kp[, 3L],
               col = "orange", size = 8)
      rgl_call("axes3d")
      rgl_call("title3d",
               xlab = "x", ylab = "y", zlab = "z",
               main = sprintf("Catmull-Rom (%s)", x$curve_type))
    }))
  }

  # ---- base-R 2D projection fallback --------------------------------------
  op <- graphics::par(mfrow = c(1L, 3L), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  proj <- function(xi, yi, xlab, ylab) {
    graphics::plot(smooth[, xi], smooth[, yi],
                   type = "l", col = col, lwd = 2,
                   xlab = xlab, ylab = ylab,
                   main = sprintf("%s-%s projection", xlab, ylab), ...)
    graphics::points(kp[, xi], kp[, yi],
                     col = "orange", pch = pch, cex = cex)
  }

  proj(1L, 2L, "x", "y")
  proj(1L, 3L, "x", "z")
  proj(2L, 3L, "y", "z")

  invisible(x)
}
