%
%  Copyright (C) 2005-2009 Friedrich Leisch
%  $Id: kcca.Rd 4650 2010-12-08 11:59:21Z leisch $
%
\name{kcca}
\alias{kcca}
\alias{kccaFamily}
\alias{flexclust-class}
\alias{kcca-class}
\alias{kccasimple-class}
\alias{kccaFamily-class}
\alias{show,kccasimple-method}
\alias{summary,kccasimple-method}
\title{K-Centroids Cluster Analysis}
\description{
  Perform k-centroids clustering on a data matrix.
}
\usage{
kcca(x, k, family=kccaFamily("kmeans"), weights=NULL, group=NULL,
     control=NULL, simple=FALSE, save.data=FALSE)
kccaFamily(which=NULL, dist=NULL, cent=NULL, name=which,
           preproc = NULL, trim=0, groupFun = "minSumClusters")

\S4method{summary}{kccasimple}(object)
}
\arguments{
  \item{x}{A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with all
    numeric columns).}
  
  \item{k}{Either the number of clusters, or a vector of cluster
    assignments, or a matrix of initial
    (distinct) cluster centroids.  If a number, a random set of (distinct)
    rows in \code{x} is chosen as the initial centroids.}
  
  \item{family}{Object of class \code{kccaFamily}.}

  \item{weights}{An optional vector of weights to be used in the clustering
    process, cannot be combined with all families.}

  \item{group}{An optional grouping vector for the data, see details
    below.}

  \item{control}{An object of class \code{flexclustControl}.}

  \item{simple}{Return an object of class \code{kccasimple}?}

  \item{save.data}{Save a copy of \code{x} in the return object?}

  \item{which}{One of \code{"kmeans"}, \code{"kmedians"},
    \code{"angle"}, \code{"jaccard"}, or \code{"ejaccard"}.}

  \item{name}{Optional long name for family, used only for show methods.}

  \item{dist}{A function for distance computation, ignored
    if \code{which} is specified.}

  \item{cent}{A function for centroid computation, ignored
    if \code{which} is specified.}

  \item{preproc}{Function for data preprocessing.}

  \item{trim}{A number in between 0 and 0.5, if non-zero then trimmed
    means are used for the \code{kmeans} family, ignored by all other
    families.}

  \item{groupFun}{Function or name of function to obtain clusters for
    grouped data, see details below.}

  \item{object}{Object of class \code{"kcca"}.}
}
\details{
  See the paper \emph{A Toolbox for K-Centroids Cluster Analysis}
  referenced below for details.
}
\section{Predefined Families}{
  Function \code{kccaFamily()} currently has the following predefined
  families (distance / centroid):
  \describe{
    \item{kmeans:}{Euclidean distance / mean}
    \item{kmedians:}{Manhattan distance / median}
    \item{angle:}{angle between observation and centroid / standardized
      mean}
    \item{jaccard:}{Jaccard distance / numeric optimization}
    \item{ejaccard:}{Jaccard distance / mean}
  }
  See Leisch (2006) for details on all combinations.
}  
\section{Group Constraints}{  
  If \code{group} is not \code{NULL}, then observations from the same
  group are restricted to belong to the same cluster (must-link
  constraint) or different clusters (cannot-link constraint) during the
  fitting process. If \code{groupFun = "minSumClusters"}, then all group
  members are 
  assign to the cluster where the center has minimal average distance to
  the group members. If \code{groupFun = "majorityClusters"}, then all
  group members are assigned to the cluster the majority would belong to
  without a constraint. 

  \code{groupFun = "differentClusters"} implements a cannot-link
  constraint, i.e., members of one group are not allowed to belong to
  the same cluster. The optimal allocation for each group is found by
  solving a linear sum assignment problem using
  \code{\link[clue]{solve_LSAP}}. Obviously the group sizes must be smaller
  than the number of clusters in this case.

  Ties are broken at random in all cases.
  Note that at the moment not all methods for fitted
  \code{"kcca"} objects respect the grouping information, most
  importantly the plot method when a data argument is specified. 
}
\references{
  Friedrich Leisch. A Toolbox for K-Centroids Cluster Analysis.
  Computational Statistics and Data Analysis, 51 (2), 526--544, 2006.

  Friedrich Leisch and Bettina Gruen. Extending standard cluster
  algorithms to allow for group constraints. In Alfredo Rizzi and
  Maurizio Vichi, editors, Compstat 2006-Proceedings in Computational
  Statistics, pages 885-892. Physica Verlag, Heidelberg, Germany, 2006.

}
\value{
  Function \code{kcca} returns objects of class \code{"kcca"} or
  \code{"kccasimple"} depending on the value of argument
  \code{simple}. The simpler objects contain fewer slots and hence are
  faster to compute, but contain no auxiliary information used by the
  plotting methods. Most plot methods for \code{"kccasimple"} objects do
  nothing and return a warning. If only centroids, cluster membership or
  prediction for new data are of interest, then the simple objects are
  sufficient.
}
\author{Friedrich Leisch}
\seealso{\code{\link{stepFlexclust}}, \code{\link{cclust}}}
\examples{
data("Nclus")
plot(Nclus)

## try kmeans 
cl1 = kcca(Nclus, k=4)
cl1

image(cl1)
points(Nclus)

## A barplot of the centroids 
barplot(cl1)


## now use k-medians and kmeans++ initialization, cluster centroids
## should be similar...

cl2 = kcca(Nclus, k=4, family=kccaFamily("kmedians"),
           control=list(initcent="kmeanspp"))
cl2

## ... but the boundaries of the partitions have a different shape
image(cl2)
points(Nclus)
}
\keyword{cluster}
