#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: conversion.R 222 2017-03-03 16:29:43Z leisch $
#

as.kcca <- function(object, ...) UseMethod("as.kcca")

as.kcca.kmeans <- function(object, data, save.data=FALSE, ...)
{
    call <- match.call()
    call[[1]] <- as.name("as.kcca")
    
    data <- as.matrix(data)
    fam <- kccaFamily("kmeans")

    z <- newKccaObject(x=data,
                                   family=fam, centers=object$centers)
    z@converged <- TRUE
    z@iter <- as.integer(1)
    z@call <- call

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=data)

    z
}

###**********************************************************

as.kcca.skmeans <- function(object, data, save.data = FALSE, ...) 
{
    call <- match.call()
    call[[1]] <- as.name("as.kcca")
    fam <- kccaFamily("angle")
    data <- fam@preproc(as.matrix(data))
    z <- newKccaObject(x = data, family = fam,
                                   centers = object$prototypes)
    z@converged <- TRUE
    z@iter <- as.integer(1)
    z@call <- call
    if (save.data) 
        z@data <- ModelEnvMatrix(designMatrix = data)
    z
}

###**********************************************************

as.kcca.partition <- function(object, data=NULL, save.data=FALSE, ...)
{
    call <- match.call()
    call[[1]] <- as.name("as.kcca")

    if(is.null(data)){
        if(is.null(object$data))
            stop("partition object does not contain data")
        else
            data <- object$data
    }
    else{
        data <- as.matrix(data)
    }

    fam <- kccaFamily("kmeans")
    if("metric" %in% names(object$call)){
        metric <- match.arg(object$call[["metric"]],
                            c("euclidean", "manhattan"))
        if(metric=="manhattan")
            fam <- kccaFamily("kmedians")
    }

    z <- newKccaObject(x=data,
                                   family=fam, centers=object$medoids)
    z@converged <- TRUE
    z@iter <- as.integer(1)
    z@call <- call

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=data)

    z
}

###**********************************************************

Cutree <- function(tree, k=NULL, h=NULL)
{
    cl <- cutree(tree, k, h)
    o <- tree$order
    ncl <- integer(length(cl))
    a <- 1:max(cl)
    aa <- 1
    while(length(a)){
        b <- min(which(ncl[o]==0))
        c <- cl[o[b]]
        ncl[cl==c] <- aa
        a <- a[a!=c]
        aa <- aa+1
    }
    ncl
}


as.kcca.hclust <- function(object, data, k, family=NULL, save.data=FALSE, ...)
{
    data <- as.matrix(data)
    call <- match.call()
    call[[1]] <- as.name("as.kcca")

    if(is.null(family) & !is.null(object$dist.method))
    {
        if(object$dist.method=="euclidean")
            family <- kccaFamily("kmeans")
        else if(object$dist.method=="manhattan")
            family <- kccaFamily("kmedians")
        else if(object$dist.method=="binary")
            family <- kccaFamily("jaccard")
    }
    
    if(is.null(family))
        stop("Cannot automatically detect correct family, please pass family argument.\n")

    cluster <- Cutree(object, k=k[1])
    centers <- family@allcent(data, cluster)
       
    z <- newKccaObject(x=data, family=family,
                                   centers=centers)

    nok <- sum(cluster != clusters(z))
    if(nok>0)
        warning(paste(nok, "cluster memberships have changed."))
    
    z@converged <- TRUE
    z@iter <- as.integer(1)
    z@call <- call

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=data)

    z
}

###**********************************************************



###**********************************************************

setAs("kccasimple", "kmeans",
function(from, to, strict=TRUE){
  if(from@family@name != "kmeans")
    warning("Family name is not kmeans.")
  
  z <- list(cluster=clusters(from),
            centers=from@centers,
            size=info(from, "size"))
  z$withinss <- double(from@k)
  for(kk in 1:from@k){
    z$withinss[kk] <- sum(from@cldist[clusters(from)==kk,1]^2)
  }
  class(z) <- "kmeans"
  z
})
  
          
  
