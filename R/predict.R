#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: predict.R 1849 2005-10-10 06:15:57Z leisch $
#

setGeneric("cluster",
           function(object, ...) standardGeneric("cluster"))

setMethod("cluster", signature(object="flexclust"),
function(object)
{
    object@cluster
})


setMethod("predict", signature(object="kcca"),
function(object, newdata=NULL, ...)
{
    if(is.null(newdata))
        return(object@cluster)
    else{
        newdata <- as(newdata, "matrix")
        newdata <- object@family@preproc(newdata)
        return(object@family@cluster(newdata, object@centers))
    }
})

