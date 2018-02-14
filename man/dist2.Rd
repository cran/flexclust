%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: dist2.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{dist2}
\alias{dist2}
\title{Compute pairwise distances between two data sets}
\description{
  This function computes and returns the distance matrix computed by
  using the specified distance measure to compute the pairwise distances
  between the rows of two data matrices.
}
\usage{
dist2(x, y, method = "euclidean", p=2)
}
\arguments{
  \item{x}{A data matrix.}
  \item{y}{A vector or second data matrix.}
  \item{method}{the distance measure to be used. This must be one of
          \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
	  \code{"canberra"}, \code{"binary"} or \code{"minkowski"}. Any
	  unambiguous substring can be given.}
  \item{p}{The power of the Minkowski distance.}
}
\details{
  This is a two-data-set equivalent of the standard function
  \code{\link[stats]{dist}}. It returns a matrix of all pairwise
  distances between rows in \code{x} and \code{y}. The current
  implementation is efficient only if \code{y} has not too many
  rows (the code is vectorized in \code{x} but not in \code{y}).
}
\note{
  The definition of Canberra distance was wrong for negative data
  prior to version 1.3-5.
}
\seealso{\code{\link[stats]{dist}}}
\author{Friedrich Leisch}
\keyword{multivariate}
\keyword{cluster}
\examples{
x = matrix(rnorm(20), ncol=4)
rownames(x) = paste("X", 1:nrow(x), sep=".")
y = matrix(rnorm(12), ncol=4)
rownames(y) = paste("Y", 1:nrow(y), sep=".")

dist2(x, y)
dist2(x, y, "man")

data(milk)
dist2(milk[1:5,], milk[4:6,])
}
