%
%  Copyright (C) 2009 Friedrich Leisch
%  $Id: randIndex.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{randIndex}
\alias{comPart}
\alias{comPart,flexclust,flexclust-method}
\alias{comPart,numeric,numeric-method}
\alias{comPart,flexclust,numeric-method}
\alias{comPart,numeric,flexclust-method}
\alias{randIndex}
\alias{randIndex,table,missing-method}
\alias{randIndex,ANY,ANY-method}
\title{Compare Partitions}
\description{
  Compute the (adjusted) Rand, Jaccard and Fowlkes-Mallows
  index for agreement of two partitions.
}
\usage{
comPart(x, y, type=c("ARI","RI","J","FM"))
\S4method{comPart}{flexclust,flexclust}(x, y, type)
\S4method{comPart}{numeric,numeric}(x, y, type)
\S4method{comPart}{flexclust,numeric}(x, y, type)
\S4method{comPart}{numeric,flexclust}(x, y, type)

randIndex(x, y, correct=TRUE, original=!correct)
\S4method{randIndex}{table,missing}(x, y, correct=TRUE, original=!correct)
\S4method{randIndex}{ANY,ANY}(x, y, correct=TRUE, original=!correct)
}
\arguments{
  \item{x}{Either a 2-dimensional cross-tabulation of cluster
    assignments (for \code{randIndex} only), an object inheriting from
    class \code{"flexclust"}, 
    or an integer vector of cluster memberships.}
  \item{y}{An object inheriting from class
    \code{"flexclust"}, or an integer vector of cluster memberships.}
  \item{type}{character vector of abbreviations of indices to compute.}
  \item{correct, original}{Logical, correct the Rand index for agreement by chance?}
}
\section{Rand Index}{
  Let \eqn{A} denote the number of all pairs of data
  points which are either put into the same cluster by both partitions or
  put into different clusters by both partitions. Conversely, let \eqn{D}
  denote the number of all pairs of data points that are put into one
  cluster in one partition, but into different clusters by the other
  partition.  The partitions disagree for all pairs \eqn{D} and
  agree for all pairs \eqn{A}. We can measure the agreement by the Rand
  index \eqn{A/(A+D)} which is invariant with respect to permutations of
  cluster labels.

  The index has to be corrected for agreement by chance if the sizes
  of the clusters are not uniform (which is usually the case), or if there
  are many clusters, see Hubert \& Arabie (1985) for details.
}
\section{Jaccard Index}{
  If the number of clusters is very large, then usually the vast
  majority of pairs of points will not be in the same cluster. The
  Jaccard index tries to account for this by using only pairs of points
  that are in the same cluster in the defintion of \eqn{A}.
}
\section{Fowlkes-Mallows}{Let \eqn{A} again be the pairs of points that
  are in the same cluster in both partitions. Fowlkes-Mallows divides
  this number by the geometric mean of the sums of the number of pairs in each
  cluster of the two partitions.  This gives the probability that a pair
  of points which are in the same cluster in one partition are also in the
  same cluster in the other partition.
}
\value{A vector of indices.}
\references{
Lawrence Hubert and Phipps Arabie. Comparing partitions.
Journal of Classification, 2, 193--218, 1985.

Marina Meila. Comparing clusterings - an axiomatic view. In Stefan
Wrobel and Luc De Raedt, editors, Proceedings of the International
Machine Learning Conference (ICML). ACM Press, 2005.
}
\author{Friedrich Leisch}
\examples{
## no class correlations: corrected Rand almost zero
g1 <- sample(1:5, size=1000, replace=TRUE)
g2 <- sample(1:5, size=1000, replace=TRUE)
tab <- table(g1, g2)
randIndex(tab)

## uncorrected version will be large, because there are many points
## which are assigned to different clusters in both cases
randIndex(tab, correct=FALSE)
comPart(g1, g2)

## let pairs (g1=1,g2=1) and (g1=3,g2=3) agree better
k <- sample(1:1000, size=200)
g1[k] <- 1
g2[k] <- 1
k <- sample(1:1000, size=200)
g1[k] <- 3
g2[k] <- 3
tab <- table(g1, g2)

## the index should be larger than before
randIndex(tab, correct=TRUE, original=TRUE)
comPart(g1, g2)
}
\keyword{cluster}
