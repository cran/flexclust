%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: info.Rd 388 2025-02-26 14:53:31Z gruen $
%
\name{info}
\alias{info,flexclust,character-method}
\title{Get Information on Fitted Flexclust Objects}
\description{
  Returns descriptive information about fitted flexclust objects like
  cluster sizes or sum of within-cluster distances.
}
\usage{
\S4method{info}{flexclust,character}(object, which, drop=TRUE, ...)
}
\arguments{
  \item{object}{Fitted object.}

  \item{which}{Which information to get. Use \code{which="help"} to list
    available information.}

  \item{drop}{Logical. If \code{TRUE} the result is coerced to the lowest
    possible dimension.}
  
  \item{\dots}{Passed to methods.}
}
\details{
  Function \code{info} can be used to access slots of fitted flexclust
  objects in a portable way, and in addition computes some
  meta-information like sum of within-cluster distances.

  Function \code{\link[modeltools]{infoCheck}} returns a logical value that is \code{TRUE}
  if the requested information can be computed from the \code{object}.
}
\seealso{\code{\link[modeltools]{info}}}
\author{Friedrich Leisch}
\examples{
data("Nclus")
plot(Nclus)

cl1 <- cclust(Nclus, k=4)
summary(cl1)

## these two are the same
info(cl1)
info(cl1, "help")

## cluster sizes
i1 <- info(cl1, "size")
i1

## average within cluster distances
i2 <- info(cl1, "av_dist")
i2

## the sum of all within-cluster distances
i3 <- info(cl1, "distsum")
i3

## sum(i1*i2) must of course be the same as i3
stopifnot(all.equal(sum(i1*i2), i3))



## This should return TRUE
modeltools::infoCheck(cl1, "size")
## and this FALSE
modeltools::infoCheck(cl1, "Homer Simpson")
## both combined
i4 <- modeltools::infoCheck(cl1, c("size", "Homer Simpson"))
i4

stopifnot(all.equal(i4, c(TRUE, FALSE)))
}
\keyword{cluster}
