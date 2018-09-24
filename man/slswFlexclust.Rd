\name{slswFlexclust}
\alias{slswFlexclust}
\alias{plot,resampleFlexclust,missing-method}
\alias{boxplot,resampleFlexclust-method}
\alias{densityplot,resampleFlexclust-method}
\alias{show,resampleFlexclust-method}
\alias{summary,resampleFlexclust-method}
\title{
    Segment Level Stability Within Solution.
}
\description{
  Assess segment level stability within solution.
}
\usage{
slswFlexclust(x, object, ...)
\S4method{plot}{resampleFlexclust,missing}(x, y, ...)
\S4method{boxplot}{resampleFlexclust}(x, which=1, ylab=NULL, ...)
\S4method{densityplot}{resampleFlexclust}(x, data, which=1, ...)
\S4method{summary}{resampleFlexclust}(object)
}
\arguments{
  \item{x}{A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with all
    numeric columns) passed to \code{\link{stepFlexclust}}.}
  \item{object}{Object of class \code{"kcca"} for \code{slwsFlexclust}
    and \code{"resampleFlexclust"} for the summary method.}
  \item{y}{Missing.}
  \item{which}{Integer or character indicating which validation measure
    is used for plotting.}
  \item{ylab}{Axis label.}
  \item{data}{Not used.}
  \item{\dots}{Additional arguments; for details see below.}
}
\details{

  Additional arguments in \code{slswFlexclust} are argument \code{nsamp}
  which is by default equal to 100 and allows to change the number of
  bootstrap pairs drawn.  Argument \code{seed} allows to set a random
  seed and argument \code{multicore} is by default \code{TRUE} and
  indicates if bootstrap samples should be drawn in parallel. Argument
  \code{verbose} is by default equal to \code{FALSE} and if \code{TRUE}
  progress information is shown during computations.

  There are plotting as well as printing and summary methods implemented
  for objects of class \code{"resampleFlexclust"}. In addition to a
  standard \code{plot} method also methods for \code{densityplot} and
  \code{boxplot} are provided.

  For more details see Dolnicar and Leisch (2017) and Dolnicar et
  al. (2018).
}
\value{
  An object of class \code{"resampleFlexclust"}.
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
\seealso{\code{\link{slsaplot}}}
\examples{
data("Nclus")
cl3 <- kcca(Nclus, k = 3)
slsw.cl3 <- slswFlexclust(Nclus, cl3, nsamp = 20)
plot(Nclus, col = clusters(cl3))
plot(slsw.cl3)
densityplot(slsw.cl3)
boxplot(slsw.cl3)
}
\keyword{cluster}
\keyword{hplot}
