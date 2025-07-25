% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometry-plane.R
\name{project_plane}
\alias{project_plane}
\title{Project plane to a surface}
\usage{
project_plane(
  target,
  width,
  height,
  shape,
  initial_positions,
  translate_first = TRUE,
  diagnostic = FALSE
)
}
\arguments{
\item{target}{target surface to be projected to, must be object that can
be converted to \code{'mesh3d'} (\code{'rgl'} surface), for example,
\code{'fs.surface'} (from \code{'freesurferformat'} package) or
\code{'ieegio_surface'} from \code{'ieegio'} package.}

\item{width, height}{width and height of the plane in world space (for
\code{'ECoG'} grid, the unit is millimeter)}

\item{shape}{vector of two integers: the first element is the number of
vertices (or electrode contacts) along \code{'width'} direction; the second
element is the number of vertices along \code{'height'} direction. The total
number of vertices of the plane will be \code{prod(shape)}. Notice}

\item{initial_positions}{a \code{shape[[1]] x shape[[2]] x 3} array or a
\code{n x 3} matrix, where \code{n} is \code{prod(shape)}, the number of
vertices indicating the initial vertex positions of the plane}

\item{translate_first}{whether to translate the plane first if the
plane center is far from the surface; default is \code{FALSE}; set to
\code{TRUE} for a warm start}

\item{diagnostic}{whether to plot diagnostic figures showing the morphing
progress.}
}
\value{
The projected vertex locations, same order as \code{initial_positions}.
}
\description{
Project a two-dimensional plane (such as \code{'ECoG'} grid) to a
three-dimensional surface while preserving the order
}
\examples{


# Construct target surface

sphere <- vcg_sphere()
target <- structure(
  class = "mesh3d",
  list(
    vb = cbind(
      sphere$vb[1:3, ] - c(0.8, 0, 0),
      sphere$vb[1:3, ] + c(0.8, 0, 0)
    ),
    it = cbind(
      sphere$it[1:3, ],
      sphere$it[1:3, ] + ncol(sphere$vb)
    )
  )
)
n_surfverts <- ncol(target$vb)

plane <- plane_geometry(width = 3, height = 3, shape = c(30, 30))
plane$vb <- plane$vb[1:3, , drop = FALSE] + c(0, 0, 2)
n_contacts <- ncol(plane$vb)

# First plot
x <- t(cbind(target$vb, plane$vb))
colnames(x) <- c('x', 'y', 'z')
graphics::pairs(
  x = x, asp = 1,
  col = c(
    rep("black", n_surfverts),
    rep("green", n_contacts)
  ),
  pch = c(
    rep(46, n_surfverts),
    rep(20, n_contacts)
  )
)

projected <- project_plane(
  target = target, width = 3, height = 3, shape = c(30, 30),
  initial_positions = t(plane$vb),
  translate_first = TRUE, diagnostic = FALSE
)

y <- rbind(x, projected)
graphics::pairs(
  x = y, asp = 1,
  col = c(
    rep("black", ncol(target$vb)),
    rep("green", n_contacts),
    rep("red", n_contacts)
  ),
  pch = c(
    rep(46, n_surfverts),
    rep(1, n_contacts),
    rep(20, n_contacts)
  )
)


}
