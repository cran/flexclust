%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: qtclust.Rd 3887 2008-02-18 10:19:40Z leisch $
%
\name{qtclust}
\alias{qtclust}
\title{Stochastic QT Clustering}
\description{
  Perform stochastic QT clustering on a data matrix.
}
\usage{
qtclust(x, radius, family = kccaFamily("kmeans"), control = NULL,
        simple=TRUE, save.data=FALSE)
}
\arguments{
  \item{x}{A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with all
    numeric columns).}

  \item{radius}{Maximum radius of clusters.}

  \item{family}{Object of class \code{kccaFamily}.}

  \item{control}{An object of class \code{flexclustControl}.}

  \item{simple}{Return an object of class \code{kccasimple}?}

  \item{save.data}{Save a copy of \code{x} in the return object?}
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
  %% A copy of this paragraph is in kcca.Rd.
  Function \code{qtclust} returns objects of class \code{"kcca"} or
  \code{"kccasimple"} depending on the value of argument
  \code{simple}. The simpler objects contain fewer slots and hence are
  faster to compute, but contain no auxiliary information used by the
  plotting methods. Most plot methods for \code{"kccasimple"} objects do
  nothing and return a warning. If only centroids, cluster membership or
  prediction for new data are of interest, then the simple objects are
  sufficient.
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
