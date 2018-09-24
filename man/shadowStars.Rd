%
%  Copyright (C) 2009 Friedrich Leisch
%  $Id: shadowStars.Rd 252 2018-09-17 08:40:24Z gruen $
%
\name{shadowStars}
\alias{shadowStars}
\alias{panelShadowViolin}
\alias{panelShadowBP}
\alias{panelShadowSkeleton}
\alias{panelShadowStripes}
\title{Shadow Stars}
\description{
  Shadow star plots and corresponding panel functions.
}
\usage{
shadowStars(object, which=1:2, project=NULL,
            width=1, varwidth=FALSE,
            panel=panelShadowStripes,
            box=NULL, col=NULL, add=FALSE, ...)

panelShadowStripes(x, col, ...)
panelShadowViolin(x, ...)
panelShadowBP(x, ...)
panelShadowSkeleton(x, ...)
}
\arguments{
  \item{object}{An object of class \code{"kcca"}.}
  \item{which}{Index numbers of dimensions of (projected) input space
    to plot.}
  \item{project}{Projection object for which a \code{predict} method
    exists, e.g., the result of \code{\link{prcomp}}.}
  \item{width}{Width of vertices connecting the cluster centroids.}
  \item{varwidth}{Logical, shall all vertices have the same width or
    should the width be proportional to number of points shown on the
    vertex?}
  \item{panel}{Function used to draw vertices.}
  \item{box}{Color of rectangle drawn around each vertex.}
  \item{col}{A vector of colors for the clusters.}
  \item{add}{Logical, start a new plot?}
  \item{...}{Passed on to panel function.}
  \item{x}{Shadow values of data points corresponding to the vertex.}
}
\details{
  
The shadow value of each data point is defined as twice the distance to
the closest centroid divided by the sum of distances to closest and
second-closest centroid. If the shadow values of a point is close to 0,
then the point is close to its cluster centroid. If the shadow value is
close to 1, it is almost equidistant to the two centroids. Thus, a
cluster that is well separated from all other clusters should have many
points with small shadow values.

The neighborhood graph of a cluster solution connects two centroids by a
vertex if at least one data point has the two centroids as closest and
second closest. The width of the vertex is proportional to the sum of
shadow values of all points having these two as closest and second
closest. A shadow star depicts the distribution of shadow values on the
vertex, see Leisch (2009) for details.

Currently four panel functions are available:
\describe{
  \item{\code{panelShadowStripes}:}{line segment for each shadow value.}
  \item{\code{panelShadowViolin}:}{violin plot of shadow values.}
  \item{\code{panelShadowBP}:}{box-percentile plot of shadow values.}
  \item{\code{panelShadowSkeleton}:}{average shadow value.}
}
}
\seealso{\code{\link{shadow}}}
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

shadowStars(c5)
shadowStars(c5, varwidth=TRUE)

shadowStars(c5, panel=panelShadowViolin)
shadowStars(c5, panel=panelShadowBP)

## always use varwidth=TRUE with panelShadowSkeleton, otherwise a few
## large shadow values can lead to misleading results:
shadowStars(c5, panel=panelShadowSkeleton)
shadowStars(c5, panel=panelShadowSkeleton, varwidth=TRUE)
}
