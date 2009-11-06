%
%  Copyright (C) 2009 Friedrich Leisch
%  $Id: randIndex.Rd 4382 2009-07-30 09:22:51Z leisch $
%
\name{randIndex}
\alias{randIndex}
\alias{randIndex,table,missing-method}
\alias{randIndex,flexclust,flexclust-method}
\alias{randIndex,integer,integer-method}
\alias{randIndex,flexclust,integer-method}
\alias{randIndex,integer,flexclust-method}
\title{Rand Index}
\description{
  Compute the Rand Index for agreement of two partitions.
}
\usage{
randIndex(x, y, correct=TRUE)
\S4method{randIndex}{table,missing}(x, y, correct=TRUE)
\S4method{randIndex}{flexclust,flexclust}(x, y, correct=TRUE)
\S4method{randIndex}{integer,integer}(x, y, correct=TRUE)
\S4method{randIndex}{flexclust,integer}(x, y, correct=TRUE)
\S4method{randIndex}{integer,flexclust}(x, y, correct=TRUE)
}
\arguments{
  \item{x}{Either a 2-dimensional cross-tabulation of cluster
    assignments, or an object inheriting from class \code{"flexclust"},
    or an integer vector of cluster memberships.}
  \item{y}{An (optional) object inheriting from class
    \code{"flexclust"}, or an integer vector of cluster memberships.}
  \item{correct}{Logical, correct the index for agreement by chance?}
  
}
\details{
Suppose we want to compare two partitions summarized by the
contingency table \eqn{T=[t_{ij}]} where \eqn{i,j=1,\ldots,K} and \eqn{t_{ij}}
denotes the number of data points which are in cluster \eqn{i} in the
first partition and in cluster \eqn{j} in the second partition. 
Let \eqn{A} denote the number of all pairs of data
points which are either put into the same cluster by both partitions or
put into different clusters by both partitions. Conversely, let \eqn{D}
denote the number of all pairs of data points that are put into one
cluster in one partition, but into different clusters by the other
partition.  The partitions disagree for all pairs \eqn{D} and
agree for all pairs \eqn{A}. We can measure the agreement by the Rand
index \eqn{A/(A+D)} which is invariant with respect to permutations of
the columns or rows of \eqn{T}.

The index has to be corrected for agreement by chance if the sizes
of the clusters are not uniform (which is usually the case), or if there
are many clusters, see Hubert \& Arabie (1985) for details.
}
\value{
  A number between -1 and 1 for the corrected version, a number between
  0 and 1 for the original version.
}
\references{
Lawrence Hubert and Phipps Arabie. Comparing partitions.
Journal of Classification, 2, 193--218, 1985.
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


## let pairs (g1=1,g2=1) and (g1=3,g2=3) agree better
k <- sample(1:1000, size=200)
g1[k] <- 1
g2[k] <- 1
k <- sample(1:1000, size=200)
g1[k] <- 3
g2[k] <- 3
tab <- table(g1, g2)

## the index should be larger than before
randIndex(tab)
randIndex(tab, correct=FALSE)
}
\keyword{cluster}
