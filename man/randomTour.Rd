%
%  Copyright (C) 2008 Friedrich Leisch
%  $Id: randomTour.Rd 3977 2008-05-09 09:43:48Z leisch $
%
\name{randomTour}
\docType{methods}
\alias{randomTour}
\alias{randomTourMatrix}
\alias{randomTour,ANY-method}
\alias{randomTour,flexclust-method}
\alias{randomTour,matrix-method}
\title{Plot a Random Tour}
\description{
  Create a series of projection plots corresponding to a random tour
  through the data. 
}
\usage{
randomTour(object, \dots)

\S4method{randomTour}{ANY}(object, \dots)
\S4method{randomTour}{matrix}(object, \dots)
\S4method{randomTour}{flexclust}(object, data=NULL, col=NULL, \dots)

randomTourMatrix(x, directions=10,
                 steps=100, sec=4, sleep = sec/steps,
                 axiscol=2, axislab=colnames(x),
                 center=NULL, radius=1, minradius=0.01, asp=1,
                 ...)
}
\arguments{
  \item{object, x}{A matrix or an object of class \code{"flexclust"}.}
  \item{data}{Data to include in plot.}
  \item{col}{Plotting colors for data points.}
  \item{directions}{Integer value, how many different directions are
    toured.}
  \item{steps}{Integer, number of steps in each direction.}
  \item{sec}{Numerical, lower bound for the number of seconds each
    direction takes.}
  \item{sleep}{Numerical, sleep for as many seconds after each picture
    has been plotted.}
  \item{axiscol}{If not \code{NULL}, then arrows are plotted for
    projections of the original coordinate axes in these colors.}
  \item{axislab}{Optional labels for the projected axes.}
  \item{center}{Center of the coordinate system to use in projected
    space. Default is the center of the plotting region.}
  \item{radius}{Relative size of the arrows.}
  \item{minradius}{Minimum radius of arrows to include.}
  \item{asp, \dots}{Passed on to \code{randomTourMatrix} and from there to \code{plot}.}
}
\details{
  Two random locations are chosen, and data then projected onto
  hyperplanes which are orthogonal to \code{step} vectors interpolating
  the two locations. The first two coordinates of the projected data are
  plotted. If \code{directions} is larger than one, then after the first
  \code{steps} plots one more random location is chosen, and the
  procedure is repeated from the current position to the
  new location, etc..

  The whole procedure is similar to a grand tour, but no attempt is made
  to optimize subsequent directions, \code{randomTour} simply chooses a random
  direction in each iteration. Use \code{rggobi} for the real thing.

  Obviously the function needs a reasonably fast computer and graphics
  device to give a smooth impression, for \code{\link{x11}} it may be
  necessary to use \code{type="Xlib"} rather than cairo.
}
\author{Friedrich Leisch}

%%% short version of examples is in tests/randomTour.R
\examples{
if(interactive()){
  par(ask=FALSE)
  randomTour(iris[,1:4], axiscol=2:5)
  randomTour(iris[,1:4], col=as.numeric(iris$Species), axiscol=4)

  x <- matrix(runif(300), ncol=3)
  x <- rbind(x, x+1, x+2)
  cl <- cclust(x, k=3, save.data=TRUE)

  randomTour(cl, center=0, axiscol="black")

  ## now use predicted cluster membership for new data as colors
  randomTour(cl, center=0, axiscol="black",
             data=matrix(rnorm(3000, mean=1, sd=2), ncol=3))
}
} 
\keyword{methods}
\keyword{hplot}
