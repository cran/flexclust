%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: qtclust.Rd 4513 2010-02-01 18:29:27Z leisch $
%
\name{qtclust}
\alias{qtclust}
\title{Stochastic QT Clustering}
\description{
  Perform stochastic QT clustering on a data matrix.
}
\usage{
qtclust(x, radius, family = kccaFamily("kmeans"), control = NULL,
        save.data=FALSE, kcca=FALSE)
}
\arguments{
  \item{x}{A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with all
    numeric columns).}

  \item{radius}{Maximum radius of clusters.}

  \item{family}{Object of class \code{\link{kccaFamily}} specifying the
    distance measure to be used.}

  \item{control}{An object of class \code{\link{flexclustControl}}
  specifying the minimum number of observations per cluster
  (\code{min.size}), and trials per iteration (\code{ntry}, see details
  below).}.

  \item{save.data}{Save a copy of \code{x} in the return object?}

  \item{kcca}{Run \code{\link{kcca}} after the QT cluster algorithm has
    converged?}
}
\details{
  This function implements a generalization of the QT clustering algorithm by
  Heyer et al. (1999), see Scharl and Leisch (2006). The only difference
  is that in each iteration not
  all possible cluster start points are considered, but only a random
  sample of size \code{control@ntry}. In most cases the resulting
  solutions are almost
  the same at a considerable speed increase, in some cases even better
  solutions are obtained than with the original algorithm. If
  \code{control@ntry} is
  set to the size of the data set, the original algorithm as proposed by
  Heyer et al. (1999) is obtained.
}
\value{
  Function \code{qtclust} by default returns objects of class
  \code{"kccasimple"}. If argument \code{kcca} is \code{TRUE}, function
    \code{kcca()} is run afterwards (initialized on the QT cluster
    solution). Data points
    not clustered by the QT cluster algorithm are omitted from the
    \code{kcca()} iterations, but filled back into the return
    object. All plot methods defined for objects of class \code{"kcca"}
    can be used.
}
\references{
  Heyer, L. J., Kruglyak, S., Yooseph, S. (1999). Exploring expression data:
  Identification and analysis of coexpressed genes. Genome Research 9,
  1106--1115.

  Theresa Scharl and Friedrich Leisch. The stochastic QT-clust
  algorithm: evaluation of stability and variance on time-course
  microarray data. In Alfredo Rizzi and Maurizio Vichi, editors,
  Compstat 2006 -- Proceedings in Computational Statistics, pages
  1015-1022. Physica Verlag, Heidelberg, Germany, 2006.
}
\author{Friedrich Leisch}
\keyword{cluster}
\examples{
x <- matrix(10*runif(1000), ncol=2)

## maximum distrance of point to cluster center is 3
cl1 <- qtclust(x, radius=3)

## maximum distrance of point to cluster center is 1
## -> more clusters, longer runtime
cl2 <- qtclust(x, radius=1)

opar <- par(c("mfrow","mar"))
par(mfrow=c(2,1), mar=c(2.1,2.1,1,1))
plot(x, col=predict(cl1), xlab="", ylab="")
plot(x, col=predict(cl2), xlab="", ylab="")
par(opar)
}
