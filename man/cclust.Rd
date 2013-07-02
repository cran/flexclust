%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: cclust.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{cclust}
\title{Convex Clustering}
\alias{cclust}
\description{
  Perform k-means clustering, hard competitive learning or neural gas on
  a data matrix.
}
\usage{
cclust(x, k, dist = "euclidean", method = "kmeans",
       weights=NULL, control=NULL, group=NULL, simple=FALSE,
       save.data=FALSE)
}
\arguments{
  \item{x}{A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with all
    numeric columns).}
  \item{k}{Either the number of clusters, or a vector of cluster
    assignments, or a matrix of initial
    (distinct) cluster centroids.  If a number, a random set of (distinct)
    rows in \code{x} is chosen as the initial centroids.}
  \item{dist}{Distance measure, one of \code{"euclidean"} (mean square
    distance) or \code{"manhattan "} (absolute distance).}
  \item{method}{Clustering algorithm: one of \code{"kmeans"},
    \code{"hardcl"} or \code{"neuralgas"}, see details below.}
  \item{weights}{An optional vector of weights to be used in the fitting
    process. Works only in combination with hard competitive learning.}
  \item{control}{An object of class \code{cclustControl}.}
  \item{group}{Currently ignored.}
  \item{simple}{Return an object of class \code{kccasimple}?}
  \item{save.data}{Save a copy of \code{x} in the return object?}
}
\details{
  This function uses the same computational engine as the earlier
  function of the same name from package `cclust'. The main difference
  is that it returns an S4 object of class \code{"kcca"}, hence all
  available methods for \code{"kcca"} objects can be used. By default
  \code{\link{kcca}} and \code{cclust} use exactly the same algorithm,
  but \code{cclust} will usually be much faster because it uses compiled
  code.
  
  If \code{dist} is \code{"euclidean"}, the distance between the cluster
  center and the data points is the Euclidian distance (ordinary kmeans
  algorithm), and cluster means are used as centroids.
  If \code{"manhattan"}, the distance between the cluster
  center and the data points is the sum of the absolute values of the
  distances, and the column-wise cluster medians are used as centroids.

  If \code{method} is \code{"kmeans"}, the classic kmeans algorithm as
  given by MacQueen (1967) is
  used, which works by repeatedly moving all cluster
  centers to the mean of their respective Voronoi sets. If
  \code{"hardcl"},
  on-line updates are used (AKA hard competitive learning), which work by
  randomly drawing an observation from \code{x} and moving the closest
  center towards that point (e.g., Ripley 1996). If
  \code{"neuralgas"} then the  neural gas algorithm by Martinetz et al
  (1993) is used. It is similar to hard competitive learning, but in
  addition to the closest centroid also the second closest centroid is
  moved in each iteration.
}
\value{
  An object of class \code{"kcca"}.
}
\seealso{
  \code{\link{cclustControl-class}}, \code{\link{kcca}}

}
\author{Evgenia Dimitriadou and Friedrich Leisch}
\references{
  MacQueen, J. (1967).  Some methods for classification and analysis of
  multivariate observations. In \emph{Proceedings of the Fifth Berkeley
  Symposium on  Mathematical Statistics and  Probability},
  eds L. M. Le Cam \& J. Neyman, 1, pp. 281--297. Berkeley, CA:
  University of California Press.

  Martinetz T., Berkovich S., and Schulten K (1993). `Neural-Gas'
  Network for Vector Quantization and its Application to Time-Series
  Prediction. IEEE Transactions on Neural Networks, 4 (4), pp. 558--569.

  Ripley, B. D. (1996)
  \emph{Pattern Recognition and Neural Networks.} Cambridge.
}
\examples{
## a 2-dimensional example
x<-rbind(matrix(rnorm(100,sd=0.3),ncol=2),
         matrix(rnorm(100,mean=1,sd=0.3),ncol=2))
cl<-cclust(x,2)
plot(x, col=predict(cl))
points(cl@centers, pch="x", cex=2, col=3) 

## a 3-dimensional example 
x<-rbind(matrix(rnorm(150,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=2,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=4,sd=0.3),ncol=3))
cl<-cclust(x, 6, method="neuralgas", save.data=TRUE)
pairs(x, col=predict(cl))
plot(cl)
}
\keyword{cluster}
