%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: pairs.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{pairs-methods}
\docType{methods}
\alias{pairs,kcca-method}
\alias{pairs,kccasimple-method}
\title{Methods for Function pairs in Package `flexclust'}
\description{
  Plot a matrix of neighbourhood graphs.
}
\usage{
  \S4method{pairs}{kcca}(x, which=NULL, project=NULL, oma=NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{"kcca"}}
  \item{which}{Index numbers of dimensions of (projected) input space
    to plot, default is to plot all dimensions.}
  \item{project}{Projection object for which a \code{predict} method
    exists, e.g., the result of \code{\link{prcomp}}.}
  \item{oma}{Outer margin.}
  \item{\dots}{Passed to the \code{plot} method.}
}
\details{
  This works only for \code{"kcca"} objects, no method is available for
  "kccasimple" objects.
  }
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
