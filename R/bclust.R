#
#  Copyright (C) 2001-2016 Friedrich Leisch
#  $Id: bclust.R 252 2018-09-17 08:40:24Z gruen $
#

bclust <- function (x, k = 2, base.iter = 10, base.k = 20, minsize = 0,
                    dist.method = "euclidian", hclust.method = "average",
                    FUN = "cclust", verbose = TRUE, final.cclust = FALSE,
                    resample=TRUE, weights=NULL, maxcluster=base.k, ...) 
{
    MYCALL <- match.call()
    x <- as.matrix(x)
    xr <- nrow(x)
    xc <- ncol(x)

    if(is.character(FUN)) FUN <- get(FUN, mode="function")

    z <- new("bclust",
             k=as.integer(k),
             allcenters = matrix(0, ncol = xc, nrow = base.iter * base.k),
             base.iter = as.integer(base.iter),
             base.k = as.integer(base.k),
             xcent = apply(x, 2, mean),
             colnames = colnames(x),
             dist.method = dist.method,
             hclust.method = hclust.method,
             maxcluster = as.integer(maxcluster),
             family=kccaFamily("kmeans"),
             call=MYCALL)
             

    optSEM <- getOption("show.error.messages")
    if(is.null(optSEM)) optSEM <- TRUE
    on.exit(options(show.error.messages = optSEM))

    if (verbose) cat("Committee Member:\n")
    for (n in 1:base.iter) {
        if (verbose){
            cat(" ", n, sep = "")
            if((n %% 20)==0) cat("\n")
        }
        if(resample){
            x1 <- x[sample(xr, replace = TRUE, prob=weights), ]
        }
        else{
            x1 <- x
        }

        for(m in 1:20){
            options(show.error.messages = FALSE)
            tryres <- try(FUN(x1, k = base.k, ...))
            if(!inherits(tryres, "try-error")) break
        }
        options(show.error.messages = optSEM)
        if(m==20){
            cat("\n",as.character(tryres),"\n")
            stop("Could not find valid cluster solution in 20 replications\n")
        }
        
        z@allcenters[((n - 1) * base.k + 1):(n * base.k),] <-
            parameters(tryres)
        
    }
    z@allcenters <-
        z@allcenters[complete.cases(z@allcenters),,drop=FALSE]
    z@allcluster <- knn1(z@allcenters, x,
                             factor(1:nrow(z@allcenters)))

    if(minsize > 0){
        z <- pruneBclust(z, x, minsize=minsize)
    }
    
    colnames(z@allcenters) <- colnames(x)
    
    if (verbose) 
        cat("\nComputing Hierarchical Clustering\n")
    z <- hclustBclust(z, x = x, k = k, final.cclust = final.cclust)
    z
}

setMethod("parameters", "bclust",
function (object, k) 
{
    centers <- matrix(0, nrow = k, ncol = ncol(object@allcenters))
    for (m in 1:k) {
        centers[m, ] <-
            colMeans(object@allcenters[object@members[,k-1] == m, , drop = FALSE])
    }
    colnames(centers) <- colnames(object@allcenters)
    centers
})

#setMethod("clusters", signature(object="bclust",newdata="missing"),
#function (object, newdata=NULL, k, ...)
#{
#    if(is.null(newdata))
#        allcluster <- object@allcluster
#    else
#        allcluster <- class::knn1(object$allcenters, newdata,
#                                  factor(1:nrow(object$allcenters)))
#    
#    return(object@members[allcluster, k - 1])
#})

setMethod("clusters", signature(object="bclust",newdata="missing"),
function (object, newdata=NULL, k=object@k, ...)
{
    allcluster <- object@allcluster
    return(object@members[allcluster, k - 1])
})


hclustBclust <-
    function (object, x, k, dist.method = object@dist.method,
              hclust.method = object@hclust.method, final.cclust = FALSE,
              maxcluster=object@maxcluster) 
{
    d <- dist(object@allcenters, method = dist.method)
    object@hclust <- stats::hclust(d, method = hclust.method)
        
    object@members <- cutree(object@hclust, 2:maxcluster)
    object@cluster <- clusters(object, k=k)
    object@centers <- parameters(object, k=k)
    if (final.cclust) {
        cclustres <- cclust(x, k = object@centers)
        object@centers <- parameters(cclustres)
        object@cluster <- clusters(cclustres)
    }
    object
}



setMethod("plot", signature(x="bclust", y="missing"),
function (x, y, maxcluster=x@maxcluster,
          main = "", ...) 
{
    opar <- par(c("mar", "oma"))
    on.exit(par(opar))

    if(nchar(main))
        par(oma = c(0, 0, 3, 0))
    
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
    par(mar = c(0, 4, 1, 1))
    plot(x@hclust, labels = FALSE, hang = -1, main="")
    x1 <- 1:maxcluster
    x2 <- 2:maxcluster
    y <- rev(x@hclust$height)[x1]
    z <- abs(diff(y))
    par(mar = c(4, 4, 1, 2))
    plot(x1, ((y - min(y))/(max(y) - min(y))), ty = "l", xlab = "", 
         ylab = "", ylim = c(0, 1))
    lines(x2, z/sum(z), col = "grey")
    text(x2, z/sum(z), labels = as.character(x2))
    mtext(main, outer = TRUE, cex = 1.5)
    layout(1)
})


### prune centers that contain not at least minsize data points
pruneBclust <- function(object, x, minsize=1, dohclust=FALSE, ...){

    ok <- FALSE
    while(!all(ok)){
        object@allcluster <- knn1(object@allcenters, x,
                                 factor(1:nrow(object@allcenters)))
        
        ok <- table(object@allcluster) >= minsize
        object@allcenters <- object@allcenters[ok, ]
    }
    if(dohclust){
        object <- hclustBclust(object, x, k = object@k, ...)
    }
    object
}

###**********************************************************

setMethod("bwplot", "bclust",
function(x, k=x@k, xlab="", strip.labels=NULL, strip.prefix="Cluster ",
         clusters=1:k, ...)
{
    if(is.null(strip.labels)){
        strip.labels <- paste(strip.prefix, 1:k, sep="")
        SIZE <- table(x@members[,k-1])
        strip.labels <-
            paste(strip.labels, ": ", SIZE, " centers", sep="")
    }

    mat <- x@allcenters
    mat <- mat[,ncol(mat):1] ## otherwise one has to read to plot from bottom to top
    DF <- matrix2df(mat, x@members[,k-1])
    DF <- DF[DF$group %in% clusters,]
    DF$group <- factor(DF$group)
    levels(DF$group) <- strip.labels[clusters]

    bwplot(variable~value|group, data=DF,
           as.table=TRUE, xlab=xlab, ...)
})
    
###**********************************************************
