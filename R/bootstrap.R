#
#  Copyright (C) 2009 Friedrich Leisch
#  $Id: bootstrap.R 4323 2009-04-04 12:44:45Z leisch $
#

bootFlexclust <- function(x, k, nboot=100, verbose=TRUE, correct=TRUE, ...)
{
    MYCALL <- match.call()
    
    nk <- length(k)
    nx <- nrow(x)

    index1 <- matrix(integer(1), nrow=nx, ncol=nboot)
    index2 <- index1
    
    clust1 <- array(integer(1), dim=c(nx, nk, nboot))
    clust2 <- clust1
    
    cent1 <- list()
    for(l in 1:nk)
    {
        cent1[[l]] <- array(double(1), dim=c(k[l], ncol(x), nboot))
    }
    cent2 <- cent1

    rand <- matrix(double(1), nrow=nboot, ncol=nk)
    colnames(rand) <- k

    for(b in 1:nboot)
    {
        if(verbose){
            if((b %% 100) == 0)
                cat("\n")
            if((b %% 10) == 0)
                cat(b, "")
        }

        index1[,b] <- i1 <- sample(1:nx, nx, replace=TRUE)
        index2[,b] <- i2 <- sample(1:nx, nx, replace=TRUE)

        s1 <- stepFlexclust(x[i1,,drop=FALSE], k=k, verbose=FALSE,
                            simple=TRUE, ...)
        s2 <- stepFlexclust(x[i2,,drop=FALSE], k=k, verbose=FALSE,
                            simple=TRUE, ...)
        
        for(l in 1:nk)
        {
            cl1 <- getModel(s1, l)
            cl2 <- getModel(s2, l)

            clust1[,l,b] <- clusters(cl1, newdata=x)
            clust2[,l,b] <- clusters(cl2, newdata=x)

            cent1[[l]][,,b] <- cl1@centers
            cent2[[l]][,,b] <- cl2@centers

            rand[b,l] <- randIndex(table(clust1[,l,b], clust2[,l,b]),
                                   correct=correct)
        }
    }
    
    if(verbose) cat("\n")

    new("bootFlexclust", k=k, centers1=cent1, centers2=cent2,
        cluster1=clust1, cluster2=clust2, index1=index1, index2=index2,
        rand=rand, call=MYCALL)
}

       
setMethod("show", "bootFlexclust",
function(object){
    cat("An object of class", sQuote(class(object)),"\n\n")
    cat("Call:\n")
    print(object@call)
    cat("\nNumber of bootstrap pairs:", nrow(object@rand),"\n")
})

setMethod("summary", "bootFlexclust",
function(object){
    cat("Call:\n")
    print(object@call)
    cat("\nSummary of Rand Indices:\n")
    print(summary(object@rand))
})

setMethod("plot", signature("bootFlexclust","missing"),
function(x, y, ...){
    boxplot(x, ...)
})

setMethod("boxplot", "bootFlexclust",
function(x, ...){
    boxplot(as.data.frame(x@rand), ...)
})

setMethod("densityplot", "bootFlexclust",
function(x, data, ...){
    Rand <- as.vector(x@rand)
    k <- rep(colnames(x@rand), rep(nrow(x@rand), ncol(x@rand)))
    k <- factor(k, levels=colnames(x@rand))

    densityplot(~Rand|k, as.table=TRUE, to=1, ...)
})





###**********************************************************

setGeneric("randIndex", function(x, y, correct=TRUE)
           standardGeneric("randIndex"))

setMethod("randIndex", signature(x="flexclust", y="flexclust"),
function(x, y, correct=TRUE){
    randIndex(table(clusters(x), clusters(y)), correct=correct)
})

setMethod("randIndex", signature(x="table", y="missing"),
function(x, y, correct=TRUE)
{
    if(length(dim(x))!=2)
        stop("Argument x needs to be a 2-dimensional table.")
    
    n <- sum(x)
    ni <- apply(x, 1, sum)
    nj <- apply(x, 2, sum)
    n2 <- choose(n, 2)

    if(correct){
        nis2 <- sum(choose(ni[ni > 1], 2))
        njs2 <- sum(choose(nj[nj > 1], 2))
        rand <- (sum(choose(x[x > 1], 2)) -
                 (nis2 * njs2)/n2)/((nis2 + njs2)/2 - (nis2 * njs2)/n2)
    }
    else
        rand <- 1 + (sum(x^2) - (sum(ni^2) + sum(nj^2))/2)/n2

    return(rand)
})

