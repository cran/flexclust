#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: kcca.R 4359 2009-06-15 16:29:23Z leisch $
#

normWeights <- function(x) x/mean(x)

mlogit <- function(x){
    x <- exp(-x)
    x/rowSums(x)
}
    
###**********************************************************


kccaFamily <- function(which=NULL, dist=NULL,
                       cent=NULL, name=which,
                       preproc=NULL,
                       trim=0, groupFun="minSumClusters")
{
    if(is.null(which) && is.null(dist))
        stop("either ", sQuote("which")," or ", sQuote("dist"),
             " must be specified\n")

    if(is.null(name)) name <- deparse(substitute(dist))

    z <- new("kccaFamily", name=name)

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

    z@cluster <- function(x, centers, n=1, distmat=NULL){

        if(is.null(distmat))
            distmat <- z@dist(x, centers)
        
        if(n==1){
            return(max.col(-distmat))
        }
        else{
            r <- t(matrix(apply(distmat, 1,
                                rank, ties.method="random"),
                          nrow=ncol(distmat)))
            z <- list()
            for(k in 1:n)
                z[[k]] <- apply(r, 1, function(x) which(x==k))
        }
        return(z)
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

    if(is.character(groupFun))
        groupFun <- get(groupFun, mode="function")
    z@groupFun <- groupFun
 
    
    z
}

###**********************************************************

kcca <- function(x, k, family=kccaFamily("kmeans"), weights=NULL,
                 group=NULL, control=NULL, simple=FALSE, save.data=FALSE)
{
    MYCALL <- match.call()
    control <- as(control, "flexclustControl")
    x <- as(x, "matrix")
    x <- family@preproc(x)
    N <- nrow(x)
    
    if(control@classify=="auto"){
        control@classify="hard"
    }

    if(!is.null(group))
    {       
        if(length(group)>N)
            warning("group vector longer than nrow(x)")
        
        group <- rep(group, length=N)
        group <- as.factor(group)
    }
    
    if(!is.null(weights))
    {
        control@classify="weighted"

        if(!family@weighted)
            stop("Centroid function of family cannot use weights")

        if(!is.null(group))
            stop("Weighted clustering for grouped data is not implemented yet.")
        ## recycle to number of observations
        weights <- rep(weights, length=N)
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
            distmat <- family@dist(x, centers)
            cluster <- family@cluster(x, distmat=distmat)

            if(control@classify == "simann"){
                cluster <- perturbClusters(cluster, sannprob)
                sannprob <- sannprob * control@simann[2]
            }

            if(!is.null(group))
                cluster <- family@groupFun(cluster, group, distmat)
            
            centers <- family@allcent(x, cluster=cluster, k=k)

            ## NAs in centers are empty clusters
            centers <- centers[complete.cases(centers),,drop=FALSE]
            k <- nrow(centers)
            
            changes <- sum(cluster!=clustold)
            if(control@verbose && (iter%%control@verbose==0))
                printIter(iter, changes, "Changes")

            if(changes==0) break
        }
    }
    else if(control@classify=="weighted")
    {
        td <- -1
        for(iter in 1:control@iter.max)
        {
            td.old <- td
            distmat <- family@dist(x, centers)
            cluster <- family@cluster(distmat=distmat)

            td <- sum(distmat[cbind(1:N, cluster)])
            if(control@verbose && (iter%%control@verbose==0))
                printIter(iter, td, "Sum of distances")
                                
            if(abs(td-td.old)/(abs(td)+0.1) < control@tolerance) break

            ## for weight computation we need similarities
            distmat <- mlogit(distmat)

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
    
    z <- newKccaObject(x=x, family=family, centers=centers, group=group,
                       iter=iter,
                       converged=(iter<control@iter.max),
                       call=MYCALL,
                       control=control,
                       simple=simple)

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=x)

    z
}

###**********************************************************

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
            k <- as.integer(k)
            if(k<2) stop("number of clusters must be at least 2")
            ## we need to avoid duplicates here
            ux <- na.omit(unique(x))
            if(nrow(ux) < k)
                stop("k larger than number of distinct complete data points in x")
            centers <- ux[sample(1:nrow(ux), k), , drop=FALSE]
            rm(ux)
        }
    }
    list(centers=centers, cluster=cluster, k=k)
}

###**********************************************************


newKccaObject <- function(x, family, centers, group=NULL, simple=FALSE,
                          ...)
{
    distmat <- family@dist(x, centers)

    z <- newKccasimpleObject(x=x, family=family, centers=centers, 
                             group=group, distmat=distmat, ...)
    
    if(!simple){
        z <- simple2kcca(x=x, from=z, group=group, distmat=distmat)
    }
    
    z
}

newKccasimpleObject <- function(x, family, centers, group=NULL,
                                distmat=NULL, ...)
{
    if(is.null(distmat))
        distmat <- family@dist(x, centers)

    cluster <- family@cluster(distmat=distmat)
    if(!is.null(group))
        cluster <- family@groupFun(cluster, group, distmat)
    names(cluster) <- rownames(x)
    colnames(centers) <- colnames(x)
    
    size <- as.vector(table(cluster))
    cldist <- as(distmat[cbind(1:nrow(x), cluster)], "matrix")
    cluster <- as(cluster, "integer")
    
    new("kccasimple",
        k=as(nrow(centers),"integer"),
        centers=centers,
        cluster=cluster,
        family=family,
        clusinfo=clusinfo(cluster, cldist, simple=TRUE),
        cldist=cldist,
        ...)
}

simple2kcca <- function(x, from, group=NULL, distmat=NULL)
{
    if(is.null(distmat))
        distmat <- from@family@dist(x, from@centers)

    cluster <- from@family@cluster(n=2, distmat=distmat)
    if(!is.null(group))
        cluster <- from@family@groupFun(cluster, group, distmat)
        
    if(ncol(distmat)>1){
        ## at least 2 clusters
        cldist <- cbind(distmat[cbind(1:nrow(x), cluster[[1]])],
                        distmat[cbind(1:nrow(x), cluster[[2]])])
        clsim <- computeClusterSim(distmat,cluster)
    }
    else{
        ## only one cluster
        cldist <- distmat
        clsim <- as.matrix(1)
    }
            
    xcent <- from@family@cent(x)
    totaldist <- sum(from@family@dist(x, matrix(xcent,nrow=1)))

    z <- as(from, "kcca")
    z@second <- as(cluster[[2]], "integer")
    z@xcent <- xcent
    z@xrange <- apply(x,2,range)
    z@totaldist <- totaldist
    z@clsim <- clsim
    z@clusinfo <- clusinfo(cluster[[1]], cldist, simple=FALSE)
    z@cldist <- cldist

    z
}

###**********************************************************


clusinfo <- function(cluster, cldist, simple=FALSE)
### cluster: vector of cluster memberships
### cldist: matrix with 1 or 2 columns
{
    size <- as.vector(table(cluster))
    
    clusinfo <-
        data.frame(size=size,
                   av_dist = as.vector(tapply(cldist[,1], cluster, sum))/size)

    if(!simple){
        clusinfo <- cbind(clusinfo,
                          max_dist = as.vector(tapply(cldist[,1], cluster, max)),
                          separation = as.vector(tapply(cldist[,2], cluster, min)))
    }
    clusinfo
}

###**********************************************************

computeClusterSim <- function(distmat, cluster)
{
    K <- max(cluster[[1]])
    z <- matrix(0, ncol=K, nrow=K)

    for(k in 1:K){
        ok1 <- cluster[[1]]==k
        if(any(ok1)){
            for(n in 1:K){
                if(k!=n){
                    ok2 <- ok1 & cluster[[2]]==n
                    if(any(ok2)){
                        z[k,n] <- 2*sum(distmat[ok2,k]/
                                        (distmat[ok2,k]+distmat[ok2,n]))
                    }
                }
            }
            z[k,] <- z[k,]/sum(ok1)
        }
    }
    diag(z) <- 1
    z
}


###**********************************************************

# make random changes in cluster membership
perturbClusters <- function(cluster, prob)
{
    x <- runif(length(cluster))
    x <- (x < prob)
    if(any(x))
        cluster[x] <- sample(1:max(cluster), sum(x), replace=TRUE)

    cluster
}


###**********************************************************

setMethod("show", "kccasimple",
function(object)
{
    cat(class(object), "object of family", sQuote(object@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    if(object@k<20){
        cat("\ncluster sizes:\n")
        print(table(clusters(object), useNA="ifany"))
    }
    else{
        cat("\n", sum(!is.na(object@cluster)), " points in ",
            object@k, " clusters", sep="")
        if(any(is.na(object@cluster))){
            cat(",", sum(is.na(object@cluster)), "outliers")
        }
        cat("\nDistribution of cluster sizes:\n", sep="")
        print(summary(as.integer(table(clusters(object)))))
    }
    cat("\n")
})
        
setMethod("summary", "kccasimple",
function(object)
{
    cat(class(object), "object of family", sQuote(object@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\ncluster info:\n")
    print(object@clusinfo)
    cat("\n")
    if(!object@converged) cat("no ")
    cat("convergence after", object@iter, "iterations\n")
    if("distsum" %in% info(object))
        cat("sum of within cluster distances:", info(object, "distsum"),"\n")
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

stepFlexclust <- function(x, k, nrep=3, verbose=TRUE,
                          FUN=kcca, drop=TRUE, group=NULL,
                          simple=FALSE, save.data=FALSE, seed=NULL,
                          multicore=TRUE, ...)
{
    MYCALL <- match.call()
    
    if(!is.null(seed)) set.seed(seed)
    
    bestKcca <- function(x, k, ...)
    {
        seed <- as.list(round(1e6*runif(nrep)))

        res <- MClapply(seed, 
                        function(y){
                            set.seed(y)
                            if(verbose) cat(" *")
                            FUN(x=x, k=k, group=group, simple=TRUE,
                                save.data=FALSE,
                                ...)
                        }, multicore=multicore)

        distsum <- sapply(res, function(y) info(y, "distsum"))
        res[[which.min(distsum)]]
    }
    
    k = as.integer(k)
    if(length(k)==0)
        return(list())
    
    z = list()
    MYCALL1 <- MYCALL
    if("drop" %in% names(MYCALL))
        MYCALL1[["drop"]] <- NULL
    
    for(n in 1:length(k)){
        if(verbose) cat(k[n], ":")
        kn <- as.character(k[n])
        z[[kn]] = bestKcca(x=x, k=k[n], ...)
        MYCALL1[["k"]] <- k[n]
        z[[kn]]@call <- MYCALL1

        if(!simple){
            ## x is usually at the beginning of kcca() pre-porcessed,
            ## here we have to do it manually!
            z[[kn]] <- simple2kcca(x=z[[kn]]@family@preproc(x),
                                   from=z[[kn]], group=group)
        }
        else{
        }
        if(verbose) cat("\n")
    }

    if(save.data){
        me <- ModelEnvMatrix(designMatrix=x)
        for(n in seq_along(z))
            z[[n]]@data <- me
    }
    
    if(drop && length(k)==1){
        return(z[[1]])
    }
    else{
        z <- new("stepFlexclust", models=z, k=as(k, "integer"),
                 nrep=as(nrep, "integer"), call=MYCALL)
        if(simple){
            x <- z@models[[1]]@family@preproc(x)
            z@xcent <- z@models[[1]]@family@cent(x)
            z@totaldist <-
                sum(z@models[[1]]@family@dist(x,
                                              matrix(z@xcent,nrow=1)))
        }
        else{
            z@xcent <- z@models[[1]]@xcent
            z@totaldist <- z@models[[1]]@totaldist
        }
        return(z)
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
    
    z <- data.frame(iter = sapply(object@models, function(x) x@iter),
                    converged = sapply(object@models, function(x) x@converged),
                    distsum = sapply(object@models,
                                     function(x) info(x, "distsum")))

    z1 <- data.frame(iter = NA,
                     converged = NA,
                     distsum = object@totaldist)
    
    z <- rbind(z1, z)
    
    print(z, na.string="")
})
    
setMethod("getModel", "stepFlexclust",
function(object, which=1)
{
    object@models[[which]]
})

setMethod("[[", signature(x="stepFlexclust", i="ANY", j="missing"),
function(x, i, j) getModel(x, i))

