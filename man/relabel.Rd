\name{relabel}
\alias{relabel}
\alias{relabel,kccasimple,character-method}
\alias{relabel,kccasimple,integer-method}
\alias{relabel,kccasimple,missing-method}
\alias{relabel,stepFlexclust,character-method}
\alias{relabel,stepFlexclust,integer-method}
\alias{relabel,stepFlexclust,missing-method}
\title{
  Relabel Cluster Results.
}
\description{
  The clusters are relabelled to obtain a unique labeling.
}
\usage{
relabel(object, by, ...)
\S4method{relabel}{kccasimple,character}(object, by, which = NULL, ...)
\S4method{relabel}{kccasimple,integer}(object, by, ...)
\S4method{relabel}{kccasimple,missing}(object, by, ...)
\S4method{relabel}{stepFlexclust,integer}(object, by = "series", ...)
\S4method{relabel}{stepFlexclust,missing}(object, by, ...)
}
\arguments{
  \item{object}{An object of class \code{"kccasimple"} or
    \code{"stepFlexclust"}.}
  \item{by}{If a character vector, it needs to be one of \code{"mean"},
    \code{"median"}, \code{"variable"}, \code{"manual"},
    \code{"centers"}, \code{"shadow"}, \code{"symmshadow"} or
    \code{"series"}.
    If missing, \code{"mean"} or \code{"series"} is used depending on if
    \code{object} is of class \code{"kccasimple"} or
    \code{"stepFlexclust"}.
    If an integer vector, it needs to indicate the new ordering.}
  \item{which}{Either an integer vector indiating the ordering or a
    vector of length one indicating the variable used for ordering.}
  \item{\dots}{Currently not used.}
}
\details{

  If \code{by} is a character vector with value \code{"mean"} or
  \code{"median"}, the clusters are ordered by the mean or median values
  over all variables for each cluster. If \code{by = "manual"}
  \code{which} needs to be a vector indicating the ordering.  If
  \code{by = "variable"} \code{which} needs to be indicate the variable
  which is used to determine the ordering. If \code{by} is
  \code{"centers"}, \code{"shadow"} or \code{"symmshadow"}, cluster
  similarities are calculated using \code{clusterSim} and used to
  determine an ordering using \code{seriate} from package
  \pkg{seriation}.

  If \code{by = "series"} the relabeling is performed over a series of
  clustering to minimize the misclassification.
}
\author{Friedrich Leisch}
\seealso{\code{\link{clusterSim}}, \code{\link[seriation]{seriate}}}
\keyword{methods}
