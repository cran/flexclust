%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: image-methods.Rd 2479 2006-03-01 15:04:04Z leisch $
%
\name{image-methods}
\docType{methods}
\alias{image,kcca-method}
\alias{image,kccasimple-method}
\title{Methods for Function image in Package `flexclust'}
\description{
  Image plot of cluster segments overlaid by neighbourhood graph.
}
\usage{
\S4method{image}{kcca}(x, which = 1:2, npoints = 100,
         xlab = "", ylab = "", fastcol = TRUE, col=NULL,
         clwd=0, graph=TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{"kcca"}.}
  \item{which}{Index number of dimensions of input space to plot.}
  \item{npoints}{Number of grid points for image.}
  \item{fastcol}{If \code{TRUE}, a greedy algorithm is used for the
    background colors of the segments, which may result in neighbouring
    segments having the same color. If \code{FALSE}, neighbouring
    segments always have different colors at a speed penalty.}
  \item{col}{Vector of background colors for the segments.}
  \item{clwd}{Line width of contour lines at cluster boundaries, use
    larger values for \code{npoints} than the default to get smooth
    lines. (Warning: We really need a smarter way to draw cluster
    boundaries!)}
  \item{graph}{Logical, add a neighborhood graph to the plot?}
  \item{xlab, ylab, \ldots}{Graphical parameters.}
}
\details{
  This works only for \code{"kcca"} objects, no method is available for
  "kccasimple" objects.
  }
\seealso{
  \code{\link{kcca}}
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
