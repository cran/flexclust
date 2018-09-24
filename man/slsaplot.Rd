\name{slsaplot}
\alias{slsaplot}
\title{
  Segment Level Stability Across Solutions Plot.
}
\description{
  Create a segment level stability across solutions plot, possibly
  using an additional variable for coloring the nodes.
}
\usage{
slsaplot(object, nodecol = NULL, ...)
}
\arguments{
  \item{object}{An object returned by \code{stepFlexclust}.}
  \item{nodecol}{A numeric vector of length equal to the number of
    observations clustered in \code{object} which represents
    an additional variable where a cluster-specific mean is calculated
    and used to color the nodes.}
  \item{\dots}{Additional graphical parameters to modify the plot.}
}
\value{
  List of length equal to the number of different cluster solutions
  minus one containing numeric vectors of the entropy values used by
  default to color the nodes.
}
\details{
  For more details see Dolnicar and Leisch (2017) and Dolnicar et
  al. (2018).
}
\references{
  Dolnicar S. and Leisch F. (2017) "Using Segment Level Stability
  to Select Target Segments in Data-Driven Market Segmentation Studies"
  Marketing Letters, 28 (3), pp. 423--436.

  Dolnicar S., Gruen B., and Leisch F. (2018) Market Segmentation
  Analysis: Understanding It, Doing It, and Making It Useful. Springer
  Singapore.
}
\author{
  Friedrich Leisch
}
\seealso{\code{\link{stepFlexclust}}, \code{\link{relabel}}, \code{\link{slswFlexclust}}}
\examples{
data("Nclus")
cl25 <- stepFlexclust(Nclus, k=2:5)
slsaplot(cl25)
cl25 <- relabel(cl25)
slsaplot(cl25)
}
\keyword{cluster}
\keyword{hplot}
