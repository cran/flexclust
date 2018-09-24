%
%  Copyright (C) 2009 Friedrich Leisch
%  $Id: stripes.Rd 252 2018-09-17 08:40:24Z gruen $
%
\name{stripes}
\alias{stripes}
\title{Stripes Plot}
\description{
  Plot distance of data points to cluster centroids using stripes.
}
\usage{
stripes(object, groups=NULL, type=c("first", "second", "all"),
        beside=(type!="first"), col=NULL, gp.line=NULL, gp.bar=NULL,
        gp.bar2=NULL, number=TRUE, legend=!is.null(groups),
        ylim=NULL, ylab="distance from centroid",
        margins=c(2,5,3,2), ...)
}
\arguments{
  \item{object}{An object of class \code{"kcca"}.}
  \item{groups}{Grouping variable to color-code the stripes. By default
    cluster membership is used as \code{groups}.}
  \item{type}{Plot distance to closest, closest and second-closest or
    to all centroids?}
  \item{beside}{Logical, make different stripes for different clusters?}
  \item{col}{Vector of colors for clusters or groups.}
  \item{gp.line, gp.bar, gp.bar2}{Graphical parameters for horizontal
    lines and background rectangular areas, see
    \code{\link[grid]{gpar}}.}
  \item{number}{Logical, write cluster numbers on x-axis?}
  \item{legend}{Logical, plot a legend for the groups?}
  \item{ylim, ylab}{Graphical parameters for y-axis.}
  \item{margins}{Margin of the plot.}
  \item{\dots}{Further graphical parameters.}
}
\details{
  A simple, yet very effective plot for visualizing the distance of each
  point from its closest and second-closest cluster centroids is a
  stripes plot. For each of the k clusters we have a rectangular area,
  which we optionally vertically
  divide into k smaller rectangles (\code{beside=TRUE}). Then we draw a
  horizontal line segment for each data point marking the distance of
  the data point from the corresponding centroid.
}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. Neighborhood graphs, stripes and shadow plots for
  cluster visualization. Statistics and Computing, 20(4), 457--469, 2010.
}
\keyword{hplot}
\examples{
bw05 <- bundestag(2005)
bavaria <- bundestag(2005, state="Bayern")

set.seed(1)
c4 <- cclust(bw05, k=4, save.data=TRUE)
plot(c4)

stripes(c4)
stripes(c4, beside=TRUE)

stripes(c4, type="sec")
stripes(c4, type="sec", beside=FALSE)
stripes(c4, type="all")

stripes(c4, groups=bavaria)

## ugly, but shows how colors of all parts can be changed
library("grid")
stripes(c4, type="all",
        gp.bar=gpar(col="red", lwd=3, fill="white"),
        gp.bar2=gpar(col="green", lwd=3, fill="black"))

}
