\name{qtclust}
\alias{qtclust}
\title{QT Clustering}
\description{
  Perform QT clustering on a data matrix.
}
\usage{
qtclust(x, radius, family = kccaFamily("kmeans"), control = NULL)
}
\arguments{
  \item{x}{A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with all
    numeric columns).}

  \item{radius}{Maximum radius of clusters.}

  \item{family}{Object of class \code{kccaFamily}.}

  \item{control}{An object of class \code{flexclustControl}.}
}
\details{
  This function implements a generalization of the QT clustering algorithm by
  Heyer et al. (1999). The only difference is that in each iteration not
  all possible cluster start points are considered, but only a random
  sample of size \code{control@ntry}. In most cases the resulting
  solutions are almost
  the same at a considerable speed increase. If \code{control@ntry} is
  set to the size of the data set, the original algorithm is obtained.
}
\value{
  An object of class \code{"kcca"}.
}
\references{
  Heyer, L. J., Kruglyak, S., Yooseph, S. (1999). Exploring expression data:
  Identification and analysis of coexpressed genes. Genome Research 9,
  1106--1115.
}
\author{Friedrich Leisch}
\keyword{cluster}
\examples{
x <- matrix(10*runif(1000), ncol=2)

## maximum distrance of point to cluster center is 3
cl1 <- qtclust(x, radius=3)

## maximum distrance of point to cluster center is 1
## -> more clusters, longer runtime
cl2 <- qtclust(x, radius=1)

opar <- par(c("mfrow","mar"))
par(mfrow=c(2,1), mar=c(2.1,2.1,1,1))
plot(x, col=predict(cl1), xlab="", ylab="")
plot(x, col=predict(cl2), xlab="", ylab="")
par(opar)
}
