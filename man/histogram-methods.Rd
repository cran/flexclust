\name{histogram-methods}
\docType{methods}
\alias{histogram,kccasimple,missing-method}
\alias{histogram,kccasimple,data.frame-method}
\alias{histogram,kccasimple,matrix-method}
\title{Methods for Function histogram in Package `flexclust'}
\description{
  Plot a histogram of the similarity of each observation to each
  cluster. 
}
\usage{
\S4method{histogram}{kccasimple,missing}(x, data, xlab="", ...)
\S4method{histogram}{kccasimple,data.frame}(x, data, xlab="", ...)
\S4method{histogram}{kccasimple,matrix}(x, data, xlab="Similarity",
          power=1, ...)
}
\arguments{
  \item{x}{An object of class \code{"kccasimple"}.}
  \item{data}{If not missing, the distance and thus similarity
    between observations and cluster centers is determined for
    the new data and used for the plots. By default the values
    from the training data are used.}
  \item{xlab}{Label for the x-axis.}
  \item{power}{Numeric indicating how similarities are transformed, for
    more details see Dolnicar et al. (2018).}
  \item{\dots}{Additional arguments passed to
    \code{\link[lattice]{histogram}}.}
}
\references{
  Dolnicar S., Gruen B., and Leisch F. (2018) Market Segmentation
  Analysis: Understanding It, Doing It, and Making It Useful. Springer
  Singapore.
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}


