\name{predict-methods}
\docType{methods}
\alias{predict,kcca-method}
\title{Predict Cluster Membership}
\description{
  Return either the cluster membership of training data or predict for
  new data.
}
\usage{
\S4method{predict}{kcca}(object, newdata, \dots)
}
\arguments{
  \item{object}{Object of class inheriting from \code{"kcca"}.}
  \item{newdata}{An optional data matrix with the same number of columns
    as the cluster centers. If omitted, the fitted values are used.}
  \item{\dots}{Currently not used.}
}
\author{Friedrich Leisch}
\keyword{methods}

