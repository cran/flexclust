#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: conversion.R 4121 2008-09-22 16:59:44Z leisch $
#

setOldClass("kmeans")
setOldClass("partition")

as.kcca <- function(object, ...) UseMethod("as.kcca")

as.kcca.kmeans <- function(object, data, save.data=FALSE, ...)
{
    data <- as.matrix(data)
    call <- match.call()
    call[[1]] <- as.name("as.kcca")
    fam <- kccaFamily("kmeans")

    z <- flexclust:::newKccaObject(x=data,
                                   family=fam, centers=object$centers)
    z@converged <- TRUE
    z@iter <- as.integer(1)
    z@call <- call

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=data)

    z
}

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

    z <- flexclust:::newKccaObject(x=data,
                                   family=fam, centers=object$medoids)
    z@converged <- TRUE
    z@iter <- as.integer(1)
    z@call <- call

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=data)

    z
}


