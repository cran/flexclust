%
%  Copyright (C) 2005-2009 Friedrich Leisch
%  $Id: stepFlexclust.Rd 4663 2011-02-04 14:07:21Z leisch $
%
\name{stepFlexclust}
\alias{stepFlexclust}
\alias{stepcclust}
\alias{stepFlexclust-class}
\alias{show,stepFlexclust-method}
\alias{plot,stepFlexclust,missing-method}
\alias{getModel}
\alias{getModel,stepFlexclust-method}
\alias{[[,stepFlexclust,ANY,missing-method}
\title{Run Flexclust Algorithms Repeatedly}
\description{
  Runs clustering algorithms repeatedly for different numbers of
  clusters and returns the minimum within cluster distance solution for
  each.
}
\usage{
stepFlexclust(x, k, nrep=3, verbose=TRUE, FUN = kcca, drop=TRUE,
              group=NULL, simple=FALSE, save.data=FALSE, seed=NULL,
              multicore=TRUE, ...)

stepcclust(...)

\S4method{plot}{stepFlexclust,missing}(x, y,
  type=c("barplot", "lines"), totaldist=NULL,
  xlab=NULL, ylab=NULL, \dots)

\S4method{getModel}{stepFlexclust}(object, which=1)
}
\arguments{
  \item{x, \dots}{passed to \code{\link{kcca}} or \code{\link{cclust}}.}

  \item{k}{A vector of integers passed in turn to the \code{k} argument
    of \code{\link{kcca}}}
  
  \item{nrep}{For each value of \code{k} run \code{\link{kcca}}
    \code{nrep} times and keep only the best solution.}

  \item{FUN}{Cluster function to use, typically \code{kcca} or
    \code{\link{cclust}}.}
  
  \item{verbose}{If \code{TRUE}, show progress information during
    computations.}
  
  \item{drop}{If \code{TRUE} and \code{K} is of length 1, then a single
    cluster object is returned instead of a \code{"stepFlexclust"}
    object.}
  
  \item{group}{An optional grouping vector for the data, see
    \code{\link{kcca}} for details.}

  \item{simple}{Return an object of class \code{kccasimple}?}

  \item{save.data}{Save a copy of \code{x} in the return object?}

  \item{seed}{If not \code{NULL}, a call to \code{set.seed()} is made
    before any clustering is done.}

  \item{multicore}{If \code{TRUE}, use package \pkg{multicore} for
    parallel processing if available.}

  \item{y}{Not used.}

  \item{type}{Create a barplot or lines plot.}

  \item{totaldist}{Include value for 1-cluster solution in plot? Default
    is \code{TRUE} if \code{K} contains \code{2}, else \code{FALSE}.}

  \item{xlab, ylab}{Graphical parameters.}
  
  \item{object}{Object of class \code{"stepFlexclust"}.}

  \item{which}{Number of model to get. If character, interpreted as
    number of clusters.}
}
\details{
  \code{stepcclust} is a simple wrapper for
  \code{stepFlexclust(...,FUN=cclust)}.

  Availability of \pkg{multicore} is checked when \pkg{flexclust} is
  loaded and stored in \code{getOption("flexclust")$have_multicore}. Set
  to \code{FALSE} for debugging and more sensible error messages in case
  something goes wrong.
}
\author{Friedrich Leisch}
\examples{
data("Nclus")
plot(Nclus)

## multicore off for CRAN checks
cl1 = stepFlexclust(Nclus, k=2:7, FUN=cclust, multicore=FALSE)
cl1

plot(cl1)

# two ways to do the same:
getModel(cl1, 4)
cl1[[4]]

opar=par("mfrow")
par(mfrow=c(2,2))
for(k in 3:6){
  image(getModel(cl1, as.character(k)), data=Nclus)
  title(main=paste(k, "clusters"))
}
par(opar)
}
\keyword{cluster}
