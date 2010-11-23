%
%  Copyright (C) 2009 Friedrich Leisch
%  $Id: bootFlexclust.Rd 4381 2009-07-24 11:48:57Z leisch $
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

  \item{multicore}{If \code{TRUE}, use package \pkg{multicore} for
    parallel processing if available. In addition, it may be a
    workstation cluster object as returned by
    \code{\link[snow]{makeCluster}} from package \pkg{snow}, see
    examples below.}
  
  \item{verbose}{If \code{TRUE}, show progress information during
    computations. Will not work with \code{multicore=TRUE}.}

  \item{y, data}{Not used.}

  \item{object}{An object of class \code{"bootFlexclust"}.}

}
\details{
  Availability of \pkg{multicore} is checked
  when \pkg{flexclust} is loaded and stored in
  \code{getOption("flexclust")$have_multicore}. Set to \code{FALSE}
  for debugging and more sensible error messages in case something
  goes wrong.
}
\author{Friedrich Leisch}
\seealso{
  \code{\link{stepFlexclust}}
}
\examples{
## data uniform on unit square
x <- matrix(runif(400), ncol=2)


## run bootstrap replications on snow cluster:
library("snow")
cl <- makeCluster(rep("localhost", 2), type = "SOCK")
clusterCall(cl, function() require("flexclust"))
bcl <- bootFlexclust(x, k=2:7, FUN=cclust, multicore=cl)
stopCluster(cl)

bcl
summary(bcl)

## splitting the square into four quadrants should be the most stable
## solution (increase nboot if not)
plot(bcl)
densityplot(bcl, from=0)
}
\keyword{cluster}
