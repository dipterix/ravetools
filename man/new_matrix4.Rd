% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-Matrix4.R
\name{new_matrix4}
\alias{new_matrix4}
\alias{as_matrix4}
\title{Create a \code{Matrix4} instance for \code{'Affine'} transform}
\usage{
new_matrix4()

as_matrix4(m)
}
\arguments{
\item{m}{a matrix or a vector to be converted to the \code{Matrix4} instance;
\code{m} must be one of the followings: for matrices, the dimension must be
\code{4x4}, \code{3x4} (the last row will be \code{0 0 0 1}), or
\code{3x3} (linear transform); for vectors, the length must be
\code{16}, \code{12} (will append \code{0 0 0 1} internally),
\code{3} (translation), or \code{1} (scale).}
}
\value{
A \code{Matrix4} instance
}
\description{
Create a \code{Matrix4} instance for \code{'Affine'} transform
}
\seealso{
\code{\link{new_vector3}}, \code{\link{new_quaternion}}
}
