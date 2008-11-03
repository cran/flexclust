%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: conversion.Rd 4123 2008-09-25 11:27:34Z leisch $
%
\name{conversion}
\alias{as.kcca}
\alias{as.kcca.kmeans}
\alias{as.kcca.partition}
\title{Convert S3 Partition Objects to KCCA}
\description{
  These functions can be used to convert the results from cluster
  functions like
  \code{\link[stats]{kmeans}} or \code{\link[cluster]{pam}} to objects
  of class \code{"kcca"}.
}
\usage{
as.kcca(object, ...)

\S3method{as.kcca}{kmeans}(object, data, save.data=FALSE, ...)
\S3method{as.kcca}{partition}(object, data=NULL, save.data=FALSE, ...)
}
\arguments{
  \item{object}{fitted object.}
  \item{data}{data which were used to obtain the clustering. For
    \code{"partition"} objects
    created by functions from package cluster this is optional, if
    \code{object} contains the data.}
  \item{save.data}{Save a copy of the data in the return object?}
  \item{\dots}{currently not used.}
}
\author{Friedrich Leisch}
\keyword{cluster}
\examples{
data(Nclus)

cl1 = kmeans(Nclus, 4)
cl1
cl1a = as.kcca(cl1, Nclus)
cl1a

\testonly{
stopifnot(all.equal(cl1$cluster, clusters(cl1a)))
}

library("cluster")
cl2 = pam(Nclus, 4)
cl2
cl2a = as.kcca(cl2)
cl2a
## the same
cl2b = as.kcca(cl2, Nclus)
cl2b

\testonly{
stopifnot(all.equal(clusters(cl2a), clusters(cl2b)))
stopifnot(all.equal(parameters(cl2a), parameters(cl2b)))
}
}
