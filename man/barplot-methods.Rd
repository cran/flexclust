%
%  Copyright (C) 2005-2009 Friedrich Leisch
%  $Id: barplot-methods.Rd 4405 2009-09-14 16:49:24Z leisch $
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

\S4method{barchart}{kcca}(x, data, xlab="",
    strip.labels=NULL, strip.prefix="Cluster ",
    col=NULL, mcol="darkred", which=NULL, legend=FALSE, ...)
}
\arguments{
  \item{height, x}{An object of class \code{"kcca"}.}
  \item{bycluster}{If \code{TRUE}, then each barplot shows one
    cluster. If \code{FALSE}, then each barplot compares all cluster for
    one input variable.}
  \item{oneplot}{If \code{TRUE}, all barplots are plotted together on
    one page, else each plot is on a separate page.}
  \item{data}{If not \code{NULL}, cluster membership is predicted for
    the new data and used for the plots. By default the
    values from the training data are used. Ignored by the \code{barchart}
    method.}
  \item{FUN}{The function to be applied to each cluster for calculating
    the bar heights. Only used, if \code{data} is not \code{NULL}.}
  \item{which}{For \code{barplot} index number of clusters for the plot,
    for \code{barchart} index numbers or names of variables to plot.}
  \item{names.arg}{A vector of names to be plotted below each bar.}
  \item{main, oma, xlab, \ldots}{Graphical parameters.}
  \item{col}{Vector of colors for the clusters.}
  \item{mcol}{If not \code{NULL}, the value of
    \code{FUN} for the complete data set is plotted over each bar as a
    line segment with color \code{mcol}.}
  \item{srt}{Number between 0 and 90, rotation of the x-axis labels.}
  \item{strip.labels}{Vector of strings for the strips of the Trellis
    display.}
  \item{strip.prefix}{Prefix string for the strips of the Trellis
    display.}
  \item{legend}{If \code{TRUE}, the barchart is always plotted on the
    current graphics device and a legend is added to the bottom of the
    plot.}
}
\note{The flexclust barchart method uses a horizontal arrangements of
    the bars, and sorts them from top to bottom. Default barcharts in
    lattice are the other way round (bottom to top). See the examples
    below how this affects, e.g., manual labels for the y axis.
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
\examples{
  cl <- cclust(iris[,-5], k=3)
  barplot(cl)
  barplot(cl, bycluster=FALSE)

  ## plot the maximum instead of mean value per cluster:
  barplot(cl, bycluster=FALSE, data=iris[,-5],
          FUN=function(x) apply(x,2,max))

  ## use lattice for plotting:
  barchart(cl)
  ## automatic abbreviation of labels
  barchart(cl, scales=list(abbreviate=TRUE))

  ## Use manual labels. Note that the flexclust barchart orders bars
  ## from top to bottom (the default does it the other way round), hence
  ## we have to rev() the labels:
  LAB <- c("SL", "SW", "PL", "PW")
  barchart(cl, scales=list(y=list(labels=rev(LAB))))
}
