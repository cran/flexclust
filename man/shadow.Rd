%
%  Copyright (C) 2009 Friedrich Leisch
%  $Id: shadow.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{shadow}
\docType{methods}
\alias{shadow}
\alias{shadow,kccasimple-method}
\alias{show,shadow-method}
\alias{plot,shadow,ANY-method}
\alias{Silhouette}
\alias{Silhouette,kcca-method}
\alias{show,Silhouette-method}
\alias{plot,Silhouette,ANY-method}
\title{Cluster shadows and silhouettes}
\description{
  Compute and plot shadows and silhouettes.
}
\usage{
\S4method{shadow}{kccasimple}(object, ...)
\S4method{Silhouette}{kcca}(object, data=NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{"kcca"} or
    \code{"kccasimple"}.}
  \item{data}{data to compute silhouette values for. If the cluster
    \code{object} was created with \code{save.data=TRUE}, then these are
    used by default.}
 \item{\dots}{currently not used.}
}
\details{
The shadow value of each data point is defined as twice the distance to
the closest centroid divided by the sum of distances to closest and
second-closest centroid. If the shadow values of a point is close to 0, then the
point is close to its cluster centroid. If the shadow value is close to 1, it
is almost equidistant to the two centroids. Thus, a cluster that is
well separated from all other clusters should have many points with
small shadow values.

The silhouette value of a data point is defined as the scaled difference
between the average dissimilarity of a point to all points in its own
cluster to the smallest average dissimilarity to the points of a
different cluster. Large silhouette values indicate good separation.

The main difference between silhouette values and shadow values is that
we replace average dissimilarities to points in a cluster by
dissimilarities to point averages (=centroids).  See Leisch (2009) for
details.
}
\seealso{\code{\link[cluster]{silhouette}}}
\references{
  Friedrich Leisch. Neighborhood graphs, stripes and shadow plots for
  cluster visualization. Statistics and Computing, 2009. Accepted for
  publication on 2009-06-16.
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
\examples{
data(Nclus)
set.seed(1)
c5 <- cclust(Nclus, 5, save.data=TRUE)
c5
plot(c5)

## high shadow values indicate clusters with *bad* separation
shadow(c5)
plot(shadow(c5))

## high Silhouette values indicate clusters with *good* separation
Silhouette(c5)
plot(Silhouette(c5))
}
