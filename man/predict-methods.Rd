%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: predict-methods.Rd 3122 2006-11-05 15:10:25Z leisch $
%
\name{predict-methods}
\docType{methods}
\alias{predict,kccasimple-method}
\alias{cluster,flexclust-method}
\title{Predict Cluster Membership}
\description{
  Return either the cluster membership of training data or predict for
  new data.
}
\usage{
\S4method{predict}{kccasimple}(object, newdata, \dots)
\S4method{cluster}{flexclust}(object)
}
\arguments{
  \item{object}{Object of class inheriting from \code{"flexclust"}.}
  \item{newdata}{An optional data matrix with the same number of columns
    as the cluster centers. If omitted, the fitted values are used.}
  \item{\dots}{Currently not used.}
}
\details{
  \code{cluster} can be used on any object of class \code{"flexclust"}
  and returns the cluster meberships of the training data.

  \code{predict} can be used only on objects of class \code{"kcca"}
  (which inherit from \code{"flexclust"}). If no \code{newdata} argument
  is specified, the function is identical to \code{cluster}, if
  \code{newdata} is specified, then cluster memberships for the new data
  are predicted.
}
\author{Friedrich Leisch}
\keyword{methods}
