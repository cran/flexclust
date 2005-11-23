%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: kcca.Rd 1779 2005-08-18 09:23:15Z leisch $
%
\name{kcca}
\alias{kcca}
\alias{kccaFamily}
\alias{flexclust-class}
\alias{kcca-class}
\alias{kccaFamily-class}
\alias{show,kcca-method}
\alias{summary,kcca-method}
\title{K-Centroids Cluster Analysis}
\description{
  Perform k-centroids clustering on a data matrix.
}
\usage{
kcca(x, k, family=kccaFamily("kmeans"), weights=NULL, group=NULL, control=NULL)
kccaFamily(which=NULL, dist=NULL, cent=NULL, name=which, similarity=FALSE,
           preproc = NULL, trim=0, groupFun = "majorityClusters")

\S4method{summary}{kcca}(object)
}
\arguments{
  \item{x}{A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with all
    numeric columns).}
  
  \item{k}{Either the number of clusters or a set of initial
    (distinct) cluster centroids.  If a number, a random set of (distinct)
    rows in \code{x} is chosen as the initial centroids.}

  \item{family}{Object of class \code{kccaFamily}.}

  \item{weights}{An optional vector of weights to be used in the clustering
    process, cannot be combined with all families.}

  \item{group}{An optional grouping vector for the data, see details
    below.}

  \item{control}{An object of class \code{flexclustControl}.}

  \item{which}{One of \code{"kmeans"}, \code{"kmedians"},
    \code{"angle"}, \code{"jaccard"}, or \code{"ejaccard"}.}

  \item{name}{Optional long name for family, used only for show methods.}

  \item{dist}{A function for distance or similarity computation, ignored
    if \code{which} is specified.}

  \item{cent}{A function for centroid computation, ignored
    if \code{which} is specified.}

  \item{similarity}{Logical, if \code{TRUE} then \code{dist} is
    interpreted as a similarity instead of a distance measure, ignored
    if \code{which} is specified.}

  \item{preproc}{Function for data preprocessing.}

  \item{trim}{A number in between 0 and 0.5, if non-zero then trimmed
    means are used for the \code{kmeans} family, ignored by all other
    families.}

  \item{groupFun}{Function or name of function to obtain clusters for
    grouped data, currently an experimental feature.}

  \item{object}{Object of class \code{"kcca"}.}
}
\details{
  If \code{group} is not \code{NULL}, then observations from the same
  group are restricted to belong to the same cluster during the fitting
  process. The cluster for each group is determined by taking the
  cluster in which the majority of the group members belong to, ties are
  broken at random. Note that at the moment not all methods for fitted
  \code{"kcca"} objects respect the grouping information, most
  importantly the plot method when a data argument is specified. 
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


## now use k-medians, cluster centroids should be similar ...
cl2 = kcca(Nclus, k=4, family=kccaFamily("kmedians"))
cl2

## ... but the boundaries of the partitions have a different shape
image(cl2)
points(Nclus)
}
\keyword{cluster}
