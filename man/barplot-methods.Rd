%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: barplot-methods.Rd 1921 2005-11-23 14:24:54Z leisch $
%
\name{barplot-methods}
\docType{methods}
\alias{barplot,kcca-method}
\title{Methods for Function barplot in Package `flexclust'}
\description{
  Barplot of cluster centers or other cluster statistics.
}
\usage{
\S4method{barplot}{kcca}(height, bycluster = TRUE, oneplot = TRUE,
    data = NULL, FUN=colMeans, main = deparse(substitute(height)), 
    which = 1:height@k, names.arg = NULL,
    oma=par("oma"), col=NULL, mcol="darkred", srt=45, ...)
}
\arguments{
  \item{height}{An object of class \code{"kcca"}.}
  \item{bycluster}{If \code{TRUE}, then each barplot shows one
    cluster. If \code{FALSE}, then each barplot compares all cluster for
    one input variable.}
  \item{oneplot}{If \code{TRUE}, all barplots are plotted together on
    one page, else each plot is on a seperate page.}
  \item{data}{If not \code{NULL}, cluster membership is predicted for
    the new data and used for the plots. By default the
    values from the traing data are used.}
  \item{FUN}{The function to be applied to each cluster for calculating
    the bar heights.}
  \item{which}{Index number of clusters for the plot.}
  \item{names.arg}{A vector of names to be plotted below each bar.}
  \item{main, oma, \ldots}{Graphical parameters.}
  \item{col}{Vector of colors for the clusters.}
  \item{mcol}{If not \code{NULL}, the value of
    \code{FUN} for the complete data set is plotted over each bar as a
    line segment with color \code{mcol}.}
  \item{srt}{Number between 0 and 90, rotation of the x-axis labels.}
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
