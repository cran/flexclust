%
%  Copyright (C) 2005-2009 Friedrich Leisch
%  $Id: propBarchart.Rd 14 2013-07-02 09:56:24Z leisch $
%
\name{propBarchart}
\alias{propBarchart}
\alias{propBarchart-class}
\alias{show,propBarchart-method}
\alias{summary,propBarchart-method}
\alias{groupBWplot}
\title{Barcharts and Boxplots for Columns of a Data Matrix Split by Groups}
\description{Split a binary or numeric matrix by a grouping variable,
  run a series of tests on all variables, adjust for multiple testing
  and graphically represent results.
}
\usage{
propBarchart(x, g, alpha=0.05, correct="holm", test="prop.test",
             sort=FALSE, strip.prefix="", strip.labels=NULL,
             which=NULL, ...)

\S4method{summary}{propBarchart}(object, ...)

groupBWplot(x, g, alpha=0.05, correct="holm", col=NULL,
            shade=!is.null(shadefun), shadefun=NULL,
            strip.prefix="", strip.labels=NULL, which=NULL, ...)
}
\arguments{
  \item{x}{A binary data matrix.}
  \item{g}{A factor specifying the groups.}
  \item{alpha}{Significance level for test of differences in
    proportions.}
  \item{correct}{Correction method for multiple testing, passed to
    \code{\link{p.adjust}}.}
  \item{test}{Test to use for detecting significant differences in
    proportions.}
  \item{sort}{Logical, sort variables by total sample mean?}
  \item{strip.prefix}{Character string prepended to strips of the
    \code{\link[lattice]{barchart}} (the remainder of the strip are group
      levels and group sizes). Ignored if \code{strip.labels} is specified.}
  \item{strip.labels}{Character vector of labels to use for strips of
    \code{\link[lattice]{barchart}}.}
  \item{which}{Index numbers or names of variables to plot.}
  \item{\ldots}{Passed on to \code{\link[lattice]{barchart}}
    or \code{\link[lattice]{bwplot}}.}
  \item{object}{Return value of \code{propBarchart}.}
  \item{col}{Vector of colors for the panels.}
  \item{shade}{If \code{TRUE}, only variables with significant 
    differences in median are filled with color.}
  \item{shadefun}{A function or name of a function to compute which
    boxes are shaded, e.g. \code{"kruskalTest"} (default),
    \code{"medianInside"} or \code{"boxOverlap"}.}
}
\details{
  
  Function \code{propBarchart} splits a binary data matrix into
  subgroups, computes the percentage of ones in each column and compares
  the proportions in the groups using \code{\link{prop.test}}. The
  p-values for all variables are adjusted for multiple testing and a
  barchart of group percentages is drawn highlighting variables with
  significant differences in proportion. The \code{summary} method can
  be used to create a corresponding table for publications.

  Function \code{groupBWplot} takes a general numeric matrix, also
  splits into subgroups and uses boxes instead of bars. By default
  \code{\link{kruskal.test}} is used to compute significant differences
  in location, in addition the heuristics from
  \code{\link{bwplot,kcca-method}} can be used. Boxes of the complete sample
  are used as reference in the background.
}
\seealso{\code{\link[flexclust]{barplot-methods}},
  \code{\link[flexclust]{bwplot,kcca-method}}}
\author{Friedrich Leisch}
\keyword{hplot}
\examples{
 ## create a binary matrix from the iris data plus a random noise column
 x <- apply(iris[,-5], 2, function(z) z>median(z))
 x <- cbind(x, Noise=sample(0:1, 150, replace=TRUE))

 ## There are significant differences in all 4 original variables, Noise
 ## has most likely no significant difference (of course the difference
 ## will be significant in alpha percent of all random samples).
 p <- propBarchart(x, iris$Species)
 p
 summary(p)
 
 x <- iris[,-5]
 x <- cbind(x, Noise=rnorm(150, mean=3))
 groupBWplot(x, iris$Species)
 groupBWplot(x, iris$Species, shade=TRUE)
 groupBWplot(x, iris$Species, shadefun="medianInside")
}
