#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: conversion.R 1849 2005-10-10 06:15:57Z leisch $
#

as.kcca <- function(object, ...) UseMethod("as.kcca")

as.kcca.kmeans <- function(object, data, ...)
{
    call <- match.call()
    call[[1]] <- as.name("as.kcca")
    fam <- kccaFamily("kmeans")

    z <- flexclust:::newKccaObject(x=data,
                                   family=fam, centers=object$centers)
    z@converged <- TRUE
    z@iter <- as.integer(1)
    z@call <- call
    z
}

as.kcca.partition <- function(object, data=NULL, ...)
{
    call <- match.call()
    call[[1]] <- as.name("as.kcca")

    if(is.null(data)){
        if(is.null(object$data))
            stop("partition object does not contain data")
        else
            data <- object$data
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
    z
}

