%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: barplot-methods.Rd 2941 2006-09-08 15:10:01Z leisch $
%
\name{barplot-methods}
\docType{methods}
\alias{barplot,kcca-method}
\alias{barplot,kccasimple-method}
\alias{barchart,kcca-method}
\alias{barchart,kccasimple-method}
\title{Barplot/chart Methods in Package `flexclust'}
\description{
  Barplot of cluster centers or other cluster statistics.
}
\usage{
\S4method{barplot}{kcca}(height, bycluster = TRUE, oneplot = TRUE,
    data = NULL, FUN=colMeans, main = deparse(substitute(height)), 
    which = 1:height@k, names.arg = NULL,
    oma=par("oma"), col=NULL, mcol="darkred", srt=45, ...)

\S4method{barchart}{kcca}(x, data, xlab="", strip.prefix="Cluster ", ...)
}
\arguments{
  \item{height, x}{An object of class \code{"kcca"}.}
  \item{bycluster}{If \code{TRUE}, then each barplot shows one
    cluster. If \code{FALSE}, then each barplot compares all cluster for
    one input variable.}
  \item{oneplot}{If \code{TRUE}, all barplots are plotted together on
    one page, else each plot is on a seperate page.}
  \item{data}{If not \code{NULL}, cluster membership is predicted for
    the new data and used for the plots. By default the
    values from the traing data are used. Ignored by the \code{barchart}
    method.}
  \item{FUN}{The function to be applied to each cluster for calculating
    the bar heights.}
  \item{which}{Index number of clusters for the plot.}
  \item{names.arg}{A vector of names to be plotted below each bar.}
  \item{main, oma, xlab, \ldots}{Graphical parameters.}
  \item{col}{Vector of colors for the clusters.}
  \item{mcol}{If not \code{NULL}, the value of
    \code{FUN} for the complete data set is plotted over each bar as a
    line segment with color \code{mcol}.}
  \item{srt}{Number between 0 and 90, rotation of the x-axis labels.}
  \item{strip.prefix}{Prefix string for the strips of the Trellis display.}
}
\details{
  This works only for \code{"kcca"} objects, no method is available for
  "kccasimple" objects.
  }
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
