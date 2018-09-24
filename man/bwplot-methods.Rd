%
%  Copyright (C) 2005-2013 Friedrich Leisch
%  $Id: barplot-methods.Rd 8 2013-06-28 16:06:50Z leisch $
%
\name{bwplot-methods}
\docType{methods}
\alias{bwplot,kcca-method}
\alias{bwplot,kccasimple-method}
\alias{bwplot,bclust-method}
\title{Box-Whisker Plot Methods in Package `flexclust'}
\description{
  Seperate boxplot of variables in each cluster in comparison with
  boxplot for complete sample.
}
\usage{
\S4method{bwplot}{kcca}(x, data, xlab="",
       strip.labels=NULL, strip.prefix="Cluster ",
       col=NULL, shade=!is.null(shadefun), shadefun=NULL, byvar=FALSE, ...)
\S4method{bwplot}{bclust}(x, k=x@k, xlab="", strip.labels=NULL, 
       strip.prefix="Cluster ", clusters=1:k, ...)
}
\arguments{
  \item{x}{An object of class \code{"kcca"} or \code{"bclust"}.}
  \item{data}{If not \code{NULL}, cluster membership is predicted for
    the new data and used for the plots. By default the
    values from the training data are used.}
  \item{xlab, \ldots}{Graphical parameters.}
  \item{col}{Vector of colors for the clusters.}
  \item{strip.labels}{Vector of strings for the strips of the Trellis
    display.}
  \item{strip.prefix}{Prefix string for the strips of the Trellis
    display.}
  \item{shade}{If \code{TRUE}, only boxes with larger 
    deviation from the median or quartiles of the total
    population of the respective variables are filled with color.}
  \item{shadefun}{A function or name of a function to compute which
    boxes are shaded, e.g. \code{"medianInside"} (default)
    or \code{"boxOverlap"}.}
  \item{byvar}{If \code{TRUE}, a panel is plotted for each variable. By
    default a panel is plotted for each group.}
  \item{k}{Number of clusters.}
  \item{clusters}{Integer vector of clusters to plot.}
}
\keyword{methods}
\keyword{hplot}
\examples{
  set.seed(1)
  cl <- cclust(iris[,-5], k=3, save.data=TRUE)
  bwplot(cl)
  bwplot(cl, byvar=TRUE)

  ## fill only boxes with color which do not contain the overall median
  ## (grey dot of background box)
  bwplot(cl, shade=TRUE)

  ## fill only boxes with color which do not overlap with the box of the
  ## complete sample (grey background box)
  bwplot(cl, shadefun="boxOverlap")
}
