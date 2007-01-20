#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: predict.R 3122 2006-11-05 15:10:25Z leisch $
#

setMethod("cluster", signature(object="flexclust"),
function(object)
{
    object@cluster
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

