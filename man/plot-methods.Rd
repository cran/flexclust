%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: plot-methods.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{plot-methods}
\docType{methods}
\alias{plot,kcca,missing-method}
\alias{plot,kccasimple,missing-method}
\title{Methods for Function plot in Package `flexclust'}
\description{
  Plot the neighbourhood graph of a cluster solution together with
  projected data points.
}
\usage{
  \S4method{plot}{kcca,missing}(x, y, which=1:2, project=NULL,
         data=NULL, points=TRUE, hull=TRUE, hull.args=NULL, 
         number = TRUE, simlines=TRUE,
         lwd=1, maxlwd=8*lwd, cex=1.5, numcol=FALSE, nodes=16,
         add=FALSE, xlab="", ylab="", xlim = NULL,
         ylim = NULL, pch=NULL, col=NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{"kcca"}}
  \item{y}{not used}
  \item{which}{Index numbers of dimensions of (projected) input space
    to plot.}
  \item{project}{Projection object for which a \code{predict} method
    exists, e.g., the result of \code{\link{prcomp}}.}
  \item{data}{Data to include in plot. If the cluster object \code{x}
    was created with \code{save.data=TRUE}, then these are used by default.}
  \item{points}{Logical, shall data points be plotted (if available)?}
  \item{hull}{If \code{TRUE}, then hulls of the data are plotted (if available).
    Can either be a logical value, one of the strings \code{"convex"}
    (the default) or \code{"ellipse"}, or a function for plotting the
    hulls.}
  \item{hull.args}{A list of arguments for the hull function.}
  \item{number}{Logical, plot number labels in nodes of graph?}
  \item{numcol, cex}{Color and size of number labels in nodes of
    graph. If \code{numcol} is logical, it switches between black and
    the color of the clusters, else it is taken as a vector of colors.}
  \item{nodes}{Plotting symbol to use for nodes if no numbers are drawn.}
  \item{simlines}{Logical, plot edges of graph?}
  \item{lwd, maxlwd}{Numerical, thickness of lines.}
  \item{add}{Logical, add to existing plot?}
  \item{xlab, ylab}{Axis labels.}
  \item{xlim, ylim}{Axis range.}
  \item{pch, col, \dots}{Plotting symbols and colors for data points.}
}
\details{
  This works only for \code{"kcca"} objects, no method is available for
  "kccasimple" objects. 
}
\references{
  Friedrich Leisch. Visualizing cluster analysis and finite mixture
  models. In Chun houh Chen, Wolfgang Haerdle, and Antony Unwin, editors,
  Handbook of Data Visualization, Springer Handbooks of Computational
  Statistics. Springer Verlag, 2008.
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
