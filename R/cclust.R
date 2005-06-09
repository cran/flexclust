cclust <- function (x, k, dist = "euclidean", method = "kmeans",
                    weights=NULL, control=NULL)
{
    MYCALL <- match.call()
    control <- as(control, "cclustControl")
    x <- as(x, "matrix")
    N <- nrow(x)

    dist <- match.arg(dist, c("euclidean", "manhattan"))
    idist <- pmatch(dist, c("euclidean", "manhattan"))

    if(dist=="euclidean"){
        family <- kccaFamily("kmeans")
    }
    else{
        family <- kccaFamily("kmedians")
    }

    method <- match.arg(method, c("kmeans", "hardcl", "neuralgas"))

    if(!is.null(weights) & (method!="hardcl")){
        warning("weights can only be used with hard competitive learning")
        method <- "hardcl"
    }

    xold <- x
    perm <- sample(N)
    x <- x[perm, ]

    centers <- initCenters(x, k, family)
    cluster <- centers$cluster
    k <- centers$k
    centers <- centers$centers


    if (method == "kmeans")
    {
        z <- .C("kmeans",
                xrows = as.integer(nrow(x)),
                xcols = as.integer(ncol(x)), 
                x = as.double(x),
                ncenters = as.integer(k), 
                centers = as.double(centers),
                cluster = as.integer(cluster), 
                iter.max = as.integer(control@iter.max),
                iter = integer(1), 
                changes = integer(control@iter.max),
                clustersize = integer(k), 
                verbose = as.integer(control@verbose),
                dist = as.integer(idist-1), PACKAGE="flexclust")
    }
    else if (method == "hardcl") {
        if(is.null(weights)){
            weights <- rep(1, N)
        }
        else{
            if(any(weights<0))
                stop("weights must be positive numbers")

            weights <- rep(weights, length.out=nrow(x))/max(weights)
            weights <- weights[perm]
        }

        if(control@method=="polynomial"){
            methrate <- 0
            rate.par <- control@pol.rate
        }
        else{
            methrate <- 1
            rate.par <- control@exp.rate
        }
        
        z <- .C("hardcl",
                xrows = as.integer(nrow(x)),
                xcols = as.integer(ncol(x)), 
                x = as.double(x),
                ncenters = as.integer(k), 
                centers = as.double(centers),
                cluster = as.integer(cluster), 
                iter.max = as.integer(control@iter.max),
                iter = integer(1), 
                clustersize = integer(k),
                verbose = as.integer(control@verbose), 
                dist = as.integer(idist-1),
                methrate = as.integer(methrate), 
                par = as.double(rate.par),
                weights = as.double(weights),
                PACKAGE="flexclust")
    }
    else if (method == "neuralgas") {
        z <- .C("neuralgas",
                xrows = as.integer(nrow(x)), 
                xcols = as.integer(ncol(x)),
                x = as.double(x),
                ncenters = as.integer(k), 
                centers = as.double(centers),
                cluster = as.integer(cluster), 
                iter.max = as.integer(control@iter.max),
                iter = integer(1), 
                clustersize = integer(k),
                verbose = as.integer(control@verbose), 
                dist = as.integer(idist-1),
                par = as.double(control@ng.rate),
                PACKAGE="flexclust")
    }

    centers <- matrix(z$centers, nrow=k)
    colnames(centers) <- colnames(x)
    
    newKccaObject(x=xold, family=family, centers=centers,
                  iter=z$iter,
                  converged=(z$iter<control@iter.max),
                  call=MYCALL,
                  control=control)
}
