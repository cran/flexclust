%
%  Copyright (C) 2005-2009 Friedrich Leisch
%  $Id: propBarchart.Rd 4358 2009-06-15 15:25:19Z leisch $
%
\name{propBarchart}
\alias{propBarchart}
\alias{propBarchart-class}
\alias{show,propBarchart-method}
\title{Barchart for Proportions of Apriori Segments in Binary Data Matrix}
\description{
  This function splits a binary data matrix into subgroups, computes the
  percantage of ones in each column and compares the proportions in the
  groups using \code{\link{prop.test}}. The p-values for all variables
  are adjusted for multiple testing and a barchart of group percentages
  is drawn highlighting variables with significant differences in
  proportion.
}
\usage{
propBarchart(x, g, alpha=0.05, correct="holm",
             strip.prefix="", strip.labels=NULL, which=NULL, ...)
}
\arguments{
  \item{x}{A binary data matrix.}
  \item{g}{A factor specifying the groups.}
  \item{alpha}{Significance level for test of differences in
    proportions.}
  \item{correct}{Correction method for multiple testing, passed to
    \code{\link{p.adjust}}.}
  \item{strip.prefix}{Character string prepended to strips of the
    \code{\link[lattice]{barchart}} (the remainder of the strip are group
      levels and group sizes). Ignored if \code{strip.labels} is specified.}
  \item{strip.labels}{Character vector of labels to use for strips of
    \code{\link[lattice]{barchart}}.}
  \item{which}{Index numbers or names of variables to plot.}
  \item{\ldots}{Passed on to \code{\link[lattice]{barchart}}.}
}
\seealso{\code{\link[flexclust]{barplot-methods}}}
\author{Friedrich Leisch}
\keyword{hplot}
\examples{
 ## create a binary matrix from the iris data plus a random noise column
 x <- apply(iris[,1:4], 2, function(z) z>median(z))
 x <- cbind(x, Noise=sample(0:1, 150, replace=TRUE))

 ## There are significant differences in all 4 original variables, Noise
 ## has most likely no significant difference (of course the difference
 ## will be significant in alpha percent of all random samples).
 propBarchart(x, iris$Species)
}
