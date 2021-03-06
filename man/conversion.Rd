%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: conversion.Rd 252 2018-09-17 08:40:24Z gruen $
%
\name{conversion}
\alias{as.kcca}
\alias{as.kcca.hclust}
\alias{as.kcca.kmeans}
\alias{as.kcca.partition}
\alias{as.kcca.skmeans}
\alias{coerce,kccasimple,kmeans-method}
\alias{Cutree}
\title{Conversion Between S3 Partition Objects and KCCA}
\description{
  These functions can be used to convert the results from cluster
  functions like
  \code{\link[stats]{kmeans}} or \code{\link[cluster]{pam}} to objects
  of class \code{"kcca"} and vice versa.
}
\usage{
as.kcca(object, ...)

\S3method{as.kcca}{hclust}(object, data, k, family=NULL, save.data=FALSE, ...)
\S3method{as.kcca}{kmeans}(object, data, save.data=FALSE, ...)
\S3method{as.kcca}{partition}(object, data=NULL, save.data=FALSE, ...)
\S3method{as.kcca}{skmeans}(object, data, save.data=FALSE, ...)
\S4method{coerce}{kccasimple,kmeans}(from, to="kmeans", strict=TRUE)

Cutree(tree, k=NULL, h=NULL)
}
\arguments{
  \item{object}{Fitted object.}
  \item{data}{Data which were used to obtain the clustering. For
    \code{"partition"} objects
    created by functions from package cluster this is optional, if
    \code{object} contains the data.}
  \item{save.data}{Save a copy of the data in the return object?}
  \item{k}{Number of clusters.}
  \item{family}{Object of class \code{"kccaFamily"}, can be omitted for
    some known distances.}
  \item{\dots}{Currently not used.}
  \item{from, to, strict}{Usual arguments for \code{\link{coerce}}}
  \item{tree}{A tree as produced by \code{\link[stats]{hclust}}.}
  \item{h}{Numeric scalar or vector with heights where the tree should
    be cut.}
}
\details{
  The standard \code{\link{cutree}} function orders clusters such that
  observation one is in cluster one, the first observation (as ordered
  in the data set) not in cluster one is in cluster two,
  etc. \code{Cutree} orders clusters as shown in the dendrogram from
  left to right such that similar clusters have similar numbers. The
  latter is used when converting to \code{kcca}. 
  
  For hierarchical clustering the cluster memberships of the converted
  object can be different from the result of \code{Cutree},
  because one KCCA-iteration has to be performed in order to obtain a
  valid \code{kcca} object. In this case a warning is issued.
}
\author{Friedrich Leisch}
\keyword{cluster}
\examples{
data(Nclus)

cl1 <- kmeans(Nclus, 4)
cl1
cl1a <- as.kcca(cl1, Nclus)
cl1a
cl1b <- as(cl1a, "kmeans")

\testonly{
stopifnot(all.equal(cl1$cluster, clusters(cl1a)))
stopifnot(all.equal(cl1$cluster, cl1b$cluster))
stopifnot(all.equal(cl1$withinss, cl1b$withinss))
}

library("cluster")
cl2 <- pam(Nclus, 4)
cl2
cl2a <- as.kcca(cl2)
cl2a
## the same
cl2b <- as.kcca(cl2, Nclus)
cl2b

\testonly{
stopifnot(all.equal(clusters(cl2a), clusters(cl2b)))
stopifnot(all.equal(parameters(cl2a), parameters(cl2b)))
}

## hierarchical clustering
hc <- hclust(dist(USArrests))
plot(hc)
rect.hclust(hc, k=3)
c3 <- Cutree(hc, k=3)
k3 <- as.kcca(hc, USArrests, k=3)
barchart(k3)
table(c3, clusters(k3))
}
