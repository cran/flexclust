%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: stepFlexclust.Rd 1779 2005-08-18 09:23:15Z leisch $
%
\name{stepFlexclust}
\alias{stepFlexclust}
\alias{stepFlexclust-class}
\alias{show,stepFlexclust-method}
\alias{plot,stepFlexclust,missing-method}
\alias{getModel}
\alias{getModel,stepFlexclust-method}
\title{Run FlexClust Algorithms Repeatedly}
\description{
  Runs clustering algorithms repeatedly for different numbers of
  clusters and returns the minimum within cluster distance solution for
  each.
}
\usage{
stepFlexclust(..., K=NULL, nrep=3, verbose=TRUE, FUN = kcca, drop=TRUE)

\S4method{plot}{stepFlexclust,missing}(x, y,
  type=c("barplot", "lines"), totaldist=NULL,
  xlab=NULL, ylab=NULL, \dots)

\S4method{getModel}{stepFlexclust}(object, which=1)
}
\arguments{
  \item{\dots}{passed to \code{\link{kcca}} or \code{\link{cclust}}.}

  \item{K}{A vector of integers passed in turn to the \code{k} argument
    of \code{\link{kcca}}}
  
  \item{nrep}{For each value of \code{k} run \code{\link{kcca}}
    \code{nrep} times and keep only the best solution.}

  \item{FUN}{Cluster function to use, typically \code{kcca} or
    \code{\link{cclust}}.}
  
  \item{verbose}{If \code{TRUE}, show progress information during
    computations.}
  
  \item{drop}{If \code{TRUE} and \code{K} is of length 1, then a single
    cluster object is returned instead of a \code{"stepFlexclust"} object.}

  \item{x}{Object of class \code{"stepFlexclust"}.}

  \item{y}{Not used.}

  \item{type}{Create a barplot or lines plot.}

  \item{totaldist}{Include value for 1-cluster solution in plot? Default
    is \code{TRUE} if \code{K} contains \code{2}, else \code{FALSE}.}

  \item{xlab, ylab}{Graphical parameters.}
  
  \item{object}{Object of class \code{"stepFlexclust"}.}

  \item{which}{Number of model to get. If character, interpreted as
    number of clusters.}
}
\author{Friedrich Leisch}
\examples{
data("Nclus")
plot(Nclus)

cl1 = stepFlexclust(Nclus, K=2:7, FUN=cclust)
cl1

plot(cl1)

getModel(cl1, 4)

opar=par("mfrow")
par(mfrow=c(2,2))
for(k in 3:6){
  image(getModel(cl1, as.character(k)), data=Nclus)
  title(main=paste(k, "clusters"))
}
par(opar)
}
\keyword{cluster}
