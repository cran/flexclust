normWeights <- function(x) x/mean(x)

mlogit <- function(x){
    x <- exp(-x)
    x/rowSums(x)
}
    
###**********************************************************


kccaFamily <- function(which=NULL, dist=NULL,
                       cent=NULL, name=which,
                       similarity=FALSE, preproc=NULL,
                       trim=0)
{
    if(is.null(which) && is.null(dist))
        stop("either ", sQuote("which")," or ", sQuote("dist"),
             " must be specified\n")

    if(is.null(name)) name <- deparse(substitute(dist))

    z <- new("kccaFamily", name=name, similarity=similarity)

    if(!is.null(preproc)) z@preproc <- preproc

    if(!is.null(which)){
        which <- match.arg(which, c("kmeans", "kmedians",
                                    "angle", "jaccard",
                                    "ejaccard"))
        if(!is.null(name)) z@name <- which

        if(which == "kmeans"){
            z@dist <- distEuclidean
            if(trim==0){
                z@cent <- function(x) colMeans(x)
                z@wcent <- function(x, weights)
                    colMeans(x*normWeights(weights))
            }
            else{
                z@cent <- function(x)
                    apply(x, 2, mean, trim=trim)
                z@wcent <- function(x, weights)
                    apply(x*normWeights(weights), 2, mean, trim=trim)
            }
            z@weighted <- TRUE
        }
        else if(which == "kmedians"){
            z@dist <- distManhattan
            z@cent <- function(x) apply(x, 2, median)
        }
        else if(which == "angle"){
            z@dist <- distAngle
            z@cent <- centAngle
            z@preproc <- function(x) x/sqrt(rowSums(x^2))
        }
        else if(which == "jaccard"){
            z@dist <- distJaccard
            z@cent <- function(x) centOptim01(x, dist=distJaccard)
        }
        else if(which == "ejaccard"){
            z@dist <- distJaccard
            z@cent <- function(x) colMeans(x)
        }
    }
    else{
        if(is.character(dist))
            z@dist <- get(dist, mode="function")
        else
            z@dist <- dist

        if(is.null(cent))
            z@cent <- function(x)
                centOptim(x, dist=dist)
        else if(is.character(cent))
            z@cent <- get(cent, mode="function")
        else
            z@cent <- cent        
    }

    if(similarity){
        z@cluster <- function(x, centers, n=1, distmat=NULL){

            if(is.null(distmat))
                distmat <- z@dist(x, centers)
            
            if(n==1){
                return(apply(distmat, 1, which.max))
            }
            else{
                r <- t(apply(- distmat, 1,
                             rank, ties.method="random"))
                z <- list()
                for(k in 1:n)
                    z[[k]] <- apply(r, 1, function(x) which(x==k))
            }
            return(z)
        }
    }
    else{
        z@cluster <- function(x, centers, n=1, distmat=NULL){

            if(is.null(distmat))
                distmat <- z@dist(x, centers)
            
            if(n==1){
                return(apply(distmat, 1, which.min))
            }
            else{
                r <- t(apply(distmat, 1,
                             rank, ties.method="random"))
                z <- list()
                for(k in 1:n)
                    z[[k]] <- apply(r, 1, function(x) which(x==k))
            }
            return(z)
        }
    }

    z@allcent <- function(x, cluster, k=max(cluster, na.rm=TRUE))
    {
        centers <- matrix(NA, nrow=k, ncol=ncol(x))
        for(n in 1:k){
            if(sum(cluster==n, na.rm=TRUE)>0){
                centers[n,] <- z@cent(x[cluster==n,,drop=FALSE])
            }
        }
        centers
    }
    
    z
}

###**********************************************************

kcca <- function(x, k, family=kccaFamily("kmeans"), weights=NULL,
                 group=NULL, control=NULL)
{
    MYCALL <- match.call()
    control <- as(control, "flexclustControl")
    x <- as(x, "matrix")
    x <- family@preproc(x)
    N <- nrow(x)
    
    if(control@classify=="auto"){
        control@classify="hard"
    }

    if(!is.null(group)) group <- rep(group, length=N)
    
    if(!is.null(weights))
    {
        control@classify="weighted"

        if(!family@weighted)
            stop("Centroid function of family cannot use weights")

        if(!is.null(group))
            stop("Weighted clustering for grouped data is not implemented yet.")
        
        if(is.null(weights))
            weights <- rep(1, N)
        else
            weights <- rep(weight, length=N)
 
    }
    
    centers <- initCenters(x, k, family)
    cluster <- centers$cluster
    k <- centers$k
    centers <- centers$centers

    sannprob <- control@simann[1]
    if(control@classify %in% c("hard", "simann"))
    {
        for(iter in 1:control@iter.max){
            
            clustold <- cluster
            cluster <- family@cluster(x, centers)

            if(control@classify == "simann"){
                cluster <- perturbClusters(cluster, sannprob)
                sannprob <- sannprob * control@simann[2]
            }

            if(!is.null(group))
                cluster <- groupedClusters(cluster, group)
                    
            centers <- family@allcent(x, cluster=cluster, k=k)

            changes <- sum(cluster!=clustold)
            if(control@verbose && (iter%%control@verbose==0))
                printIter(iter, changes, "Changes")

            if(changes==0) break
        }
    }
    else if(control@classify=="weighted")
    {

        ts <- -1
        for(iter in 1:control@iter.max){

            ts.old <- ts
            distmat <- family@dist(x, centers)
            cluster <- family@cluster(distmat=distmat)

            ## for weighted kmeans we need similarities
            if(!family@similarity){
                distmat <- mlogit(distmat)
            }
            ts <- sum(distmat[cbind(1:N, cluster)])

            if(control@verbose && (iter%%control@verbose==0))
                printIter(iter, ts, "Sum of similarities")
                                
            if(abs(ts-ts.old)/(abs(ts)+0.1) < control@tolerance) break
            
            for(n in 1:k){
                w <- weights*distmat[,n]
                w[cluster==n] <- w[cluster==n]+control@gamma
                centers[n,] <- family@wcent(x, weights=w)
            }
        }
        ## make sure the final centers are centroids of their partitions
        for(n in 1:k){
            centers[n,] <- family@cent(x[cluster==n,,drop=FALSE])
        }
    }
    else
        stop("Unknown classification method")

    centers <- centers[complete.cases(centers),,drop=FALSE]
    colnames(centers) <- colnames(x)
    
    newKccaObject(x=x, family=family, centers=centers, group=group,
                  iter=iter,
                  converged=(iter<control@iter.max),
                  call=MYCALL,
                  control=control)
}


initCenters <- function(x, k, family)
{
    cluster <- integer(nrow(x))
    if(is.matrix(k)){
        centers <- k
        if(any(duplicated(centers)))
            stop("initial centers are not distinct")
        k <- nrow(k)
    }
    else{
        if(length(k)>1){
            cluster <- rep(k, length=nrow(x))
            k <- max(cluster)
            centers <- family@allcent(x, cluster=cluster, k=k)
        }
        else{
            ## we need to avoid duplicates here
            ux <- unique(x)
            if(nrow(ux) < k)
                stop("more centroids than distinct data points requested")
            centers <- ux[sample(1:nrow(ux), k), , drop=FALSE]
            rm(ux)
        }
    }
    list(centers=centers, cluster=cluster, k=k)
}


newKccaObject <- function(x, family, centers, group=NULL, ...)
{
    distmat <- family@dist(x, centers)
    cluster <- family@cluster(n=2, distmat=distmat)

    if(!is.null(group))
        cluster <- groupedClusters(cluster, group)

    if(ncol(distmat)>1){
        ## at least 2 clusters
        cldist <- cbind(distmat[cbind(1:nrow(x), cluster[[1]])],
                        distmat[cbind(1:nrow(x), cluster[[2]])])
        clsim <- clusterSim(distmat,cluster)
    }
    else{
        ## only one cluster
        cldist <- distmat
        clsim <- as.matrix(1)
    }
        
    withindist <- as.vector(tapply(cldist[,1], cluster[[1]], sum))
    
    z <- new("kcca", k=as(nrow(centers),"integer"),
             centers=centers,
             cluster=as.integer(cluster[[1]]),
             second=as.integer(cluster[[2]]),
             family=family, 
             xrange=apply(x,2,range),
             xcent=family@cent(x),
             size=table(cluster[[1]]),
             withindist=withindist,
             clsim=clsim,
             cldist=cldist, ...)
    z
}
       

clusterSim <- function(distmat, cluster)
{
    K <- max(cluster[[1]])
    z <- matrix(0, ncol=K, nrow=K)

    for(k in 1:K){
        for(n in 1:K){
            if(k!=n){
                ok <- cluster[[1]]==k & cluster[[2]]==n
                if(any(ok)){
                    z[k,n] <- 2*mean(distmat[ok,k]/
                                     (distmat[ok,k]+distmat[ok,n]))
                }
            }
        }
    }
    diag(z) <- 1
    z
}

# make random changes in cluster membership
perturbClusters <- function(cluster, prob)
{
    x <- runif(length(cluster))
    x <- (x < prob)
    if(any(x))
        cluster[x] <- sample(1:max(cluster), sum(x), replace=TRUE)

    cluster
}

## Assign each observation to the cluster where the majority from its
## group belongs to
groupedClusters <- function(cluster, group)
{
    getMax <- function(a)
    {
        b
    }

    if(is.list(cluster))
    {
        K <- max(unlist(cluster))
        ## use factors to make sure all possible clusters appear in
        ## both tables
        x <- 2*table(group, factor(cluster[[1]], levels=1:K)) +
               table(group, factor(cluster[[2]], levels=1:K))
    }
    else{
        x <- table(group, cluster)
    }

    m <- max.col(x)
    names(m) <- row.names(x)
    z <- m[group]
    names(z) <- NULL
    
    if(is.list(cluster))
    {
        ## get second best
        x[cbind(1:nrow(x), m)] <- 0
        m <- max.col(x)
        names(m) <- row.names(x)
        z1 <- m[group]
        names(z1) <- NULL
        z <- list(z, z1)
    }
    z
}
    
                


###**********************************************************

setMethod("show", "kcca",
function(object)
{
    cat("kcca object of family", sQuote(object@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\ncluster sizes:\n")
    ## the size slot does not contain the NA count for qtclust objects
    print(table(object@cluster, exclude=NULL))
    cat("\n")
    if(!object@converged) cat("no ")
    cat("convergence after", object@iter, "iterations\n")
    cat("sum of within cluster",
        ifelse(object@family@similarity, "similarities", "distances"),
        ":", sum(object@withindist), "\n")
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

###**********************************************************

setGeneric("neighbours",
function(object, symm=TRUE, ...) standardGeneric("neighbours"))

setMethod("neighbours", signature(object="kcca"),
function(object, symm=TRUE, ...)
{
    neighbours(list(object@cluster, object@second))
})

setMethod("neighbours", signature(object="list"),
function(object, symm=TRUE, ...)
{
    z <- list()
    for(n in 1:max(object[[1]], object[[2]])){
        if(symm)
            z[[n]] <- unique(c(object[[2]][object[[1]]==n],
                               object[[1]][object[[2]]==n]))
        else
            z[[n]] <- unique(object[[2]][object[[1]]==n])
    }
    z
})

###**********************************************************


setGeneric("ntable", function(object, ...)
           standardGeneric("ntable"))

setMethod("ntable", signature(object="kcca"),
function(object, freq=TRUE, distance=FALSE, ...)
{
    z <- table(object@cluster, object@second)

    if(distance){
        freq <- FALSE
        z[z==0] <- NA
        z <- length(object@cluster)-z
        diag(z) <- 0
    }
    
    if(!freq)
        z <- z/rowSums(z,na.rm=TRUE)

    z
})
        
    
###**********************************************************        

stepFlexclust <- function(..., K=NULL, nrep=3, verbose=TRUE,
                          FUN=kcca, drop=TRUE)
{

    MYCALL <- match.call()
    
    bestKcca <- function(...)
    {
        distsum <- NULL
        for(m in 1:nrep){
            if(verbose) cat(" *")
            x = FUN(...)
            if(is.null(distsum)){
                if(x@family@similarity)
                    distsum <- -Inf
                else
                    distsum <- Inf
            }

            if(x@family@similarity){
                if(sum(x@withindist) > distsum){
                    z <- x
                    distsum <- sum(x@withindist)
                }
            }
            else{
                if(sum(x@withindist) < distsum){
                    z <- x
                    distsum <- sum(x@withindist)
                }
            }
        }
        z
    }

    if(is.null(K)){
        z = bestKcca(...)
        z@call <- MYCALL
        if(verbose) cat("\n")
        return(z)
    }
    
    K = as.integer(K)
    if(length(K)==0)
        return(list())
    
    z = list()
    MYCALL1 <- MYCALL
    if("drop" %in% names(MYCALL))
        MYCALL1[["drop"]] <- NULL
    
    for(n in 1:length(K)){
        if(verbose) cat(K[n], ":")
        z[[as.character(K[n])]] = bestKcca(..., k=K[n])
        MYCALL1[["K"]] <- K[n]
        
        z[[as.character(K[n])]]@call <- MYCALL1
        if(verbose) cat("\n")
    }
    if(drop && length(K)==1){
        return(z[[1]])
    }
    else{
        return(new("stepFlexclust", models=z, K=as(K, "integer"),
                   nrep=as(nrep, "integer"), call=MYCALL))
    }
}


setMethod("show", "stepFlexclust",
function(object)
{
    cat("stepFlexclust object of family",
        sQuote(object@models[[1]]@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    
    z1 <- sapply(object@models, function(x) x@iter)
    z2 <- sapply(object@models, function(x) x@converged)

    z <- data.frame(iter=z1, converged=z2)
    if("withindist" %in% slotNames(object@models[[1]]))
        z$distsum <- sapply(object@models, function(x) sum(x@withindist))
    
    print(z)

})
    
setGeneric("getModel",
function(object, ...) standardGeneric("getModel"))

setMethod("getModel", "stepFlexclust",
function(object, which=1)
{
    object@models[[which]]
})
