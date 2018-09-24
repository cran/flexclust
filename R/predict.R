#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: predict.R 222 2017-03-03 16:29:43Z leisch $
#

setMethod("clusters", signature(object="flexclust", newdata="missing"),
function(object, newdata)
{
    object@cluster
})

setMethod("clusters", signature(object="flexclust", newdata="ANY"),
function(object, newdata, ...)
{
    predict(object, newdata, ...)
})


setMethod("predict", signature(object="kccasimple"),
function(object, newdata=NULL, ...)
{
    if(is.null(newdata)){
        z <- object@cluster
    }
    else{
        newdata <- as(newdata, "matrix")
        newdata <- object@family@preproc(newdata)
        z <- object@family@cluster(newdata, object@centers)
        names(z) <- rownames(newdata)
    }
    z
})

