%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: projAxes.Rd 1907 2005-11-15 12:08:59Z leisch $
%
\name{projAxes}
\alias{projAxes}
\alias{projAxes}
\alias{projAxes-class}
\alias{plot,projAxes,missing-method}
\alias{show,projAxes-method}
\alias{placeLabels}
\alias{placeLabels,projAxes-method}
\title{Add Arrows for Projeted Axes to a Plot}
\description{
  Adds arrows for original coordinate axes to a projection
  plot. Currently only a method for the results of
  \code{\link[stats]{prcomp}} exists.
}
\usage{
projAxes(object, which=1:2, center=NULL,
                     col="red", radius=NULL,
                     minradius=0.1, textargs=list(col=col),
                     col.names=getColnames(object),
                     which.names="", ...)

placeLabels(object)
\S4method{placeLabels}{projAxes}(object)
}
\arguments{
  \item{object}{Return value of a projection method like
    \code{\link[stats]{prcomp}}.}
  \item{which}{Index number of dimensions of (projected) input space
    that have been plotted.}
  \item{center}{Center of the coordinate system to use in projected
    space. Default is the center of the plotting region.}
  \item{col}{Color of arrows.}
  \item{radius}{Relative size of the arrows.}
  \item{minradius}{Minimum radius of arrows to include (relative to
    arrow size).}
  \item{textargs}{List of arguments for \code{\link{text}}.}
  \item{col.names}{Variable names of the original data.}
  \item{which.names}{A regular expression which variable names to
    include in the plot.}
  \item{\dots}{Passed to  \code{\link[graphics]{arrows}}.}
}
\value{
  \code{projAxes} returns an object of class \code{"projAxes"}, which is
  added to the current plot by its \code{show} and \code{plot}
  method.}
\author{Friedrich Leisch}
\keyword{hplot}
\examples{
data(milk)
milk.pca <- prcomp(milk, scale=TRUE)

## create a biplot step by step
plot(predict(milk.pca), type="n")
text(predict(milk.pca), rownames(milk), col="green", cex=0.8)
projAxes(milk.pca)

## the same, but arrows are blue, centered at origin and all arrows are
## plotted 
plot(predict(milk.pca), type="n")
text(predict(milk.pca), rownames(milk), col="green", cex=0.8)
projAxes(milk.pca, col="blue", center=0, minradius=0)

## use points instead of text, plot PC2 and PC3, manual radius
## specification, store result
plot(predict(milk.pca)[,c(2,3)])
arr <- projAxes(milk.pca, which=c(2,3), radius=1.2)
plot(arr)
}
