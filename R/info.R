#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: info.R 3122 2006-11-05 15:10:25Z leisch $
#

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
        if(all(c("size","av_dist") %in% INFOS))
            return(sum(object@clusinfo$size *
                       object@clusinfo$av_dist))
        else
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

setMethod("parameters", signature(object="kccasimple"),
function(object, ...)
{
    object@centers
})

