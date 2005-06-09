\name{plot-methods}
\docType{methods}
\alias{plot,kcca,missing-method}
\title{Methods for Function plot in Package `flexclust'}
\description{
  Plot the neighbourhood graph of a cluster solution together with
  projected data points.
}
\usage{
  \S4method{plot}{kcca,missing}(x, y, which=1:2, project=NULL,
         data=NULL, points=!is.null(data), hull=!is.null(data),
         number = TRUE, simlines=TRUE,
         lwd=1, maxlwd=8*lwd, cex=1.5, numcol="black", nodes=16,
         add=FALSE, xlab="", ylab="", xlim = NULL,
         ylim = NULL, pch=NULL, col=NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{"kcca"}}
  \item{y}{not used}
  \item{which}{Index number of dimensions of (projected) input space
    to plot.}
  \item{project}{Projection object, for which a \code{predict} method
    exists, e.g., the result of \code{\link{prcomp}}.}
  \item{data}{Data to include in plot.}
  \item{points}{Logical, shall data points be plotted?}
  \item{hull}{If \code{TRUE}, then hulls of confidence regions of
    the data are plotted. If a function, it is used for plotting the
    hulls.}
  \item{number}{Logical, plot number labels in nodes of graph?}
  \item{numcol, cex}{Color and size of number labels in nodes of graph.}
  \item{nodes}{Plotting symbol to use for nodes if no numbers are drawn.}
  \item{simlines}{Logical, plot edges of graph?}
  \item{lwd, maxlwd}{Numerical, thickness of lines.}
  \item{add}{Logical, add to existing plot?}
  \item{xlab, ylab}{Axis labels.}
  \item{xlim, ylim}{Axis range.}
  \item{pch, col, \dots}{Plotting symbols and colors for data points.}
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
