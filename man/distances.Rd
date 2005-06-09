\name{distances}
\alias{centAngle}
\alias{centOptim01}
\alias{centOptim}
\alias{distAngle}
\alias{distCor}
\alias{distEuclidean}
\alias{distJaccard}
\alias{distManhattan}
\alias{distMax}
\title{Distance and Centroid Computation}
\description{
  Helper functions to create \code{\link{kccaFamily}} objects.
}
\usage{
distAngle(x, centers)
distCor(x, centers)
distEuclidean(x, centers)
distJaccard(x, centers)
distManhattan(x, centers)
distMax(x, centers)

centAngle(x)

centOptim(x, dist)
centOptim01(x, dist)
}
\arguments{
  \item{x}{A data matrix.}
  \item{centers}{A matrix of centroids.}
  \item{dist}{A distance function.}
}
\author{Friedrich Leisch}
\keyword{cluster}
