%
%  Copyright (C) 2009 Friedrich Leisch
%  $Id: bootFlexclust.Rd 252 2018-09-17 08:40:24Z gruen $
%
\name{bootFlexclust}
\alias{bootFlexclust}
\alias{bootFlexclust-class}
\alias{show,bootFlexclust-method}
\alias{summary,bootFlexclust-method}
\alias{plot,bootFlexclust,missing-method}
\alias{boxplot,bootFlexclust-method}
\alias{densityplot,bootFlexclust-method}
\title{Bootstrap Flexclust Algorithms}
\description{
  Runs clustering algorithms repeatedly for different numbers of
  clusters on bootstrap replica of the original data and returns
  corresponding cluster assignments, centroids and Rand indices
  comparing pairs of partitions.
}
\usage{
bootFlexclust(x, k, nboot=100, correct=TRUE, seed=NULL,
              multicore=TRUE, verbose=FALSE, ...)

\S4method{summary}{bootFlexclust}(object)
\S4method{plot}{bootFlexclust,missing}(x, y, ...)
\S4method{boxplot}{bootFlexclust}(x, ...)
\S4method{densityplot}{bootFlexclust}(x, data, ...)
}
\arguments{
  \item{x, k, \ldots}{Passed to \code{\link{stepFlexclust}}.}
  
  \item{nboot}{Number of bootstrap pairs of partitions.}
  
  \item{correct}{Logical, correct the index for agreement by chance?}

  \item{seed}{If not \code{NULL}, a call to \code{set.seed()} is made
    before any clustering is done.}

  \item{multicore}{If \code{TRUE}, use package \pkg{parallel} for
    parallel processing. In addition, it may be a workstation cluster
    object as returned by \code{\link[parallel]{makeCluster}}, see
    examples below.}
  
  \item{verbose}{If \code{TRUE}, show progress information during
    computations. Will not work with \code{multicore=TRUE}.}

  \item{y, data}{Not used.}

  \item{object}{An object of class \code{"bootFlexclust"}.}

}
\details{
  Availability of \pkg{multicore} is checked
  when \pkg{flexclust} is loaded. This information is stored and can be
  obtained using
  \code{getOption("flexclust")$have_multicore}. Set to \code{FALSE}
  for debugging and more sensible error messages in case something
  goes wrong.
}
\author{Friedrich Leisch}
\seealso{
  \code{\link{stepFlexclust}}
}
\examples{
\dontrun{

## data uniform on unit square
x <- matrix(runif(400), ncol=2)

cl <- FALSE

## to run bootstrap replications on a workstation cluster do the following:
library("parallel")
cl <- makeCluster(2, type = "PSOCK")
clusterCall(cl, function() require("flexclust"))


## 50 bootstrap replicates for speed in example,
## use more for real applications
bcl <- bootFlexclust(x, k=2:7, nboot=50, FUN=cclust, multicore=cl)

bcl
summary(bcl)

## splitting the square into four quadrants should be the most stable
## solution (increase nboot if not)
plot(bcl)
densityplot(bcl, from=0)
}}
\keyword{cluster}
