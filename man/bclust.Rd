\name{bclust}
\alias{bclust}
\alias{plot,bclust,missing-method}
\alias{clusters,bclust,missing-method}
\alias{parameters,bclust-method}
\title{Bagged Clustering}
\usage{
bclust(x, k = 2, base.iter = 10, base.k = 20, minsize = 0,
       dist.method = "euclidian", hclust.method = "average",
       FUN = "cclust", verbose = TRUE, final.cclust = FALSE,
       resample = TRUE, weights = NULL, maxcluster = base.k, ...)
\S4method{plot}{bclust,missing}(x, y, maxcluster = x@maxcluster, main = "", ...)
\S4method{clusters}{bclust,missing}(object, newdata, k, ...)
\S4method{parameters}{bclust}(object, k)
}
\arguments{
  \item{x}{Matrix of inputs (or object of class \code{"bclust"} for plot).}
  \item{k}{Number of clusters.}
  \item{base.iter}{Number of runs of the base cluster algorithm.}
  \item{base.k}{Number of centers used in each repetition of the
    base method.}
  \item{minsize}{Minimum number of points in a base cluster.}
  \item{dist.method}{Distance method used for the hierarchical
    clustering, see \code{\link[stats]{dist}} for available distances.}
  \item{hclust.method}{Linkage method used for the hierarchical
    clustering, see \code{\link[stats]{hclust}} for available methods.}
  \item{FUN}{Partitioning cluster method used as base algorithm.}
  \item{verbose}{Output status messages.}
  \item{final.cclust}{If \code{TRUE}, a final cclust step is performed
    using the output of the bagged clustering as initialization.}
  \item{resample}{Logical, if \code{TRUE} the base method is run on
    bootstrap samples of \code{x}, else directly on \code{x}.}
  \item{weights}{Vector of length \code{nrow(x)}, weights for the
    resampling. By default all observations have equal weight.}
  \item{maxcluster}{Maximum number of clusters memberships are to be
    computed for.}
  \item{object}{Object of class \code{"bclust"}.}
  \item{main}{Main title of the plot.}
  \item{\dots}{Optional arguments top be passed to the base method
    in \code{bclust}, ignored in \code{plot}.}
  \item{y}{Missing.}
  \item{newdata}{An optional data matrix with the same number of columns
    as the cluster centers. If omitted, the fitted values are used.}
}
\description{
    Cluster the data in \code{x} using the bagged clustering
    algorithm. A partitioning cluster algorithm such as
    \code{\link{cclust}} is run repeatedly on bootstrap samples from the
    original data. The resulting cluster centers are then combined using
    the hierarchical cluster algorithm \code{\link[stats]{hclust}}.
}
\details{
    First, \code{base.iter} bootstrap samples of the original data in
    \code{x} are created by drawing with replacement. The base cluster
    method is run on each of these samples with \code{base.k}
    centers. The \code{base.method} must be the name of a partitioning
    cluster function returning an object with the same slots as the
    return value of \code{\link{cclust}}.

    This results in a collection of \code{iter.base * base.centers}
    centers, which are subsequently clustered using the hierarchical
    method \code{\link[stats]{hclust}}. Base centers with less than
    \code{minsize} points in there respective partitions are removed
    before the hierarchical clustering.  The resulting dendrogram is
    then cut to produce \code{k} clusters.
}
\value{
    \code{bclust} returns objects of class
    \code{"bclust"} including the slots
    \item{hclust}{Return value of the hierarchical clustering of the
	collection of base centers (Object of class \code{"hclust"}).}
    \item{cluster}{Vector with indices of the clusters the inputs are
	assigned to.}
    \item{centers}{Matrix of centers of the final clusters. Only useful,
	if the hierarchical clustering method produces convex clusters.}
    \item{allcenters}{Matrix of all \code{iter.base * base.centers}
	centers found in the base runs.}
}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. Bagged clustering. Working Paper 51, SFB ``Adaptive
  Information Systems and Modeling in Economics and Management
  Science'', August 1999. \doi{10.57938/9b129f95-b53b-44ce-a129-5b7a1168d832}

  Sara Dolnicar and Friedrich Leisch. Winter tourist segments in
  Austria: Identifying stable vacation styles using bagged clustering
  techniques. Journal of Travel Research, 41(3):281-292, 2003.
}
\seealso{\code{\link[stats]{hclust}}, \code{\link{cclust}}}
\keyword{multivariate}
\keyword{cluster}
\examples{
data(iris)
bc1 <- bclust(iris[,1:4], 3, base.k=5)
plot(bc1)

table(clusters(bc1, k=3))
parameters(bc1, k=3)
}

