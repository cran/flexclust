%
%  Copyright (C) 2008 Friedrich Leisch
%  $Id: clusterSim.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{clusterSim}
\alias{clusterSim}
\alias{clusterSim,kcca-method}
\alias{clusterSim,kccasimple-method}
\title{Cluster Similarity Matrix}
\description{
  Returns a matrix of cluster similarities. Currently two methods for
  computing similarities of clusters are implemented, see details below.
}
\usage{
\S4method{clusterSim}{kcca}(object, data=NULL, method=c("shadow", "centers"), 
           symmetric=FALSE, ...)
\S4method{clusterSim}{kccasimple}(object, data=NULL, method=c("shadow", "centers"), 
           symmetric=FALSE, ...)
}
\arguments{
  \item{object}{fitted object.}
  
  \item{data}{Data to use for computation of the shadow values. If
    the cluster object \code{x} was created with \code{save.data=TRUE},
    then these are used by default. Ignored if \code{method="centers"}.}
  
  \item{method}{Type of similarities, see details below.}

  \item{symmetric}{Compute symmetric or asymmetric shadow values?
    Ignored if \code{method="centers"}.}
  
  \item{\dots}{currently not used.}
}

\details{

  If \code{method="shadow"} (the default), then the similarity of two
  clusters is proportional to the number of points in a cluster, where
  the centroid of the other cluster is second-closest. See Leisch (2006,
  2008) for detailed formulas.

  If \code{method="centers"}, then first the pairwise distances between
  all centroids are computed and rescaled to [0,1]. The similarity
  between tow clusters is then simply 1 minus the rescaled distance.
}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. A Toolbox for K-Centroids Cluster Analysis.
  Computational Statistics and Data Analysis, 51 (2), 526--544, 2006.

  Friedrich Leisch. Visualizing cluster analysis and finite mixture
  models. In Chun houh Chen, Wolfgang Haerdle, and Antony Unwin, editors,
  Handbook of Data Visualization, Springer Handbooks of Computational
  Statistics. Springer Verlag, 2008.
}
\examples{
example(Nclus)

clusterSim(cl)
clusterSim(cl, symmetric=TRUE)

## should have similar structure but will be numerically different:
clusterSim(cl, symmetric=TRUE, data=Nclus[sample(1:550, 200),])

## different concept of cluster similarity
clusterSim(cl, method="centers")
}
\keyword{cluster}
