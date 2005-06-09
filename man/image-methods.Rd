\name{image-methods}
\docType{methods}
\alias{image,kcca-method}
\title{Methods for Function image in Package `flexclust'}
\description{
  Image plot of cluster segments overlaid by neighbourhood graph.
}
\usage{
\S4method{image}{kcca}(x, which = 1:2, npoints = 100,
         xlab = "", ylab = "", fastcol = TRUE, col=NULL, ...)
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
  \item{xlab, ylab, \ldots}{Graphical parameters.}
}
\seealso{
  \code{\link{kcca}}
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
