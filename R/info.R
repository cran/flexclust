#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: info.R 1849 2005-10-10 06:15:57Z leisch $
#

setGeneric("info",
function(object, which, ...) standardGeneric("info"))

setMethod("info", signature(object="ANY", which="missing"),
function(object, which, ...)
{
    info(object, which="help")
})
    
setMethod("info", signature(object="flexclust", which="character"),
function(object, which, drop=TRUE, ...)
{
    INFOS <- c("distsum",names(object@clusinfo))
    if("help" %in% which){
        return(INFOS)
    }

    which <- INFOS[pmatch(which, INFOS)]
    
    if(any(is.na(which))){
        stop(paste("Requested info not available, use",
                   sQuote(paste("which=",dQuote("help"),sep="")),
                   "to list available infos."))
    }

    if(any(which=="distsum")){
        return(sum(object@cldist[,1]))
    }
    else{
        z <- object@clusinfo[,which,drop=drop]
        if(is.vector(z))
            names(z) <- rownames(object@clusinfo)
        return(z)
    }
})

###**********************************************************

infoCheck <- function(object, which, ...)
{
    which %in% info(object, "help")
}

###**********************************************************

setGeneric("parameters",
           function(object, ...) standardGeneric("parameters"))

setMethod("parameters", signature(object="kcca"),
function(object, ...)
{
    object@centers
})

