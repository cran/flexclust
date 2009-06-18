%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: predict-methods.Rd 4366 2009-06-18 13:43:44Z leisch $
%
\name{predict-methods}
\docType{methods}
\alias{predict,kccasimple-method}
\alias{clusters,flexclust,missing-method}
\alias{clusters,flexclust,ANY-method}
\title{Predict Cluster Membership}
\description{
  Return either the cluster membership of training data or predict for
  new data.
}
\usage{
\S4method{predict}{kccasimple}(object, newdata, \dots)
\S4method{clusters}{flexclust,ANY}(object, newdata, \dots)
}
\arguments{
  \item{object}{Object of class inheriting from \code{"flexclust"}.}
  \item{newdata}{An optional data matrix with the same number of columns
    as the cluster centers. If omitted, the fitted values are used.}
  \item{\dots}{Currently not used.}
}
\details{
  \code{clusters} can be used on any object of class \code{"flexclust"}
  and returns the cluster memberships of the training data.

  \code{predict} can be used only on objects of class \code{"kcca"}
  (which inherit from \code{"flexclust"}). If no \code{newdata} argument
  is specified, the function is identical to \code{clusters}, if
  \code{newdata} is specified, then cluster memberships for the new data
  are predicted. \code{clusters(object, newdata, ...)} is an alias for
  \code{predict(object, newdata, ...)}.
}
\author{Friedrich Leisch}
\keyword{methods}
