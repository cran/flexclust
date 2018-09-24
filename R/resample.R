bootTT <- function(n) list(train1=sample(1:n, replace=TRUE),
                           train2=sample(1:n, replace=TRUE),
                           test=1:n)

randEval <- function(c1,c2)
    doRandIndex(table(c1,c2), correct=TRUE, original=TRUE)

bootScheme <- new("resampleScheme",
                  traintest = bootTT,
                  validate = randEval,
                  valname = c("ARI", "RI"))


###**********************************************************




clusterwiseScheme <- function(object)
{
    stopifnot(is(object, "kcca"))
    k <- object@k
    c0 <- object@cluster

    clusterwise <- function(c1,c2)
    {
        z1 <- z2 <- matrix(double(1), nrow=k, ncol=k)
        for(m in 1:k){
            ok0 <- (c0==m)
            for(n in 1:k){
                ok1 <- (c1==n)
                ok2 <- (c2==n)
                z1[m,n] <- sum(ok0&ok1)/sum(ok0|ok1)
                z2[m,n] <- sum(ok0&ok2)/sum(ok0|ok2)
            }
        }
        z1 <- apply(z1, 1, max)
        z2 <- apply(z2, 1, max)
        (z1+z2)/2
    }
        
    new("resampleScheme",
        traintest = bootTT,
        validate = clusterwise,
        valname = as.character(1:k))
}

###**********************************************************

predstr_splitSet <- function(n){
  index <- 1:n
  train <- sort(sample(x = index, size = n/2))
  test <- index[-train]
  list(train1 = train, train2=test, test = test)
}

predstr <- function(pred, testclusters){
    
  x <- length(pred)
  l <- length(unique(pred))
  
  res <- .C(C_predstr_calcSums,
            pred = as.integer(pred),
            testclusters = as.integer(testclusters),
            x = as.integer(x),
            l = as.integer(l),
            res = integer(l))[["res"]]
 
  npoints <- table(testclusters)
  for(z in 1:l){
    n <- npoints[z]
    res[z] <- (res[z] * 1)/(n*(n-1))
  }

  sort(res)[1:2]
}

predstrScheme <- new("resampleScheme",
                     traintest = predstr_splitSet,
                     validate = predstr,
                     valname = paste("PS", 1:2))


###**********************************************************

resampleFlexclust <- function(x, k, scheme=bootScheme, nsamp=100, seed=NULL,
                              multicore=TRUE, verbose=FALSE, ...)
{
    MYCALL <- match.call()

    if(!is.null(seed)) set.seed(seed)
    seed <- round(2^31 * runif(nsamp, -1, 1))

    nk <- length(k)
    nx <- nrow(x)

    ## empirical experiments show parallization does not pay for this
    ## (sample is too fast)
    index1 <- index2 <- indextest <- list()    
    for(b in 1:nsamp){
        ind <- scheme@traintest(nx)
        index1[[b]] <- ind$train1
        index2[[b]] <- ind$train2
        indextest[[b]] <- ind$test
    }

    BFUN <- function(b){
        if(verbose){
            if((b %% 100) == 0)
                cat("\n")
            if((b %% 10) == 0)
                cat(b, "")
        }

        set.seed(seed[b])
        s1 <- stepFlexclust(x[index1[[b]],,drop=FALSE], k=k, verbose=FALSE,
                            simple=TRUE, multicore=FALSE, ...)
        s2 <- stepFlexclust(x[index2[[b]],,drop=FALSE], k=k, verbose=FALSE,
                            simple=TRUE, multicore=FALSE, ...)

        clust1 <- clust2 <- matrix(integer(1), nrow=nx, ncol=nk)
        cent1 <- cent2 <- list()
        valmat <- matrix(0, ncol=nk, nrow=length(scheme@valname))
        
        for(l in 1:nk)
        {
            if(nk>1){
                cl1 <- getModel(s1, l)
                cl2 <- getModel(s2, l)
            }
            else{
                cl1 <- s1
                cl2 <- s2
            }

            clust1[,l] <- clusters(cl1, newdata=x[indextest[[b]],])
            clust2[,l] <- clusters(cl2, newdata=x[indextest[[b]],])

            cent1[[l]] <- cl1@centers
            cent2[[l]] <- cl2@centers

            valmat[,l] <- scheme@validate(clust1[,l], clust2[,l])
        }
        list(cent1=cent1, cent2=cent2, clust1=clust1, clust2=clust2,
             valmat=valmat)
    }
    
    z <- MClapply(as.list(1:nsamp), BFUN, multicore=multicore)

    clust1 <- lapply(z, function(x) x$clust1)
    clust2 <- lapply(z, function(x) x$clust2)

    cent1 <- cent2 <- list()
    for(l in 1:nk){
        cent1[[l]] <- unlist(lapply(z, function(x) x$cent1[[l]]))
        cent2[[l]] <- unlist(lapply(z, function(x) x$cent2[[l]]))
        dim(cent1[[l]]) <- dim(cent2[[l]]) <- c(k[l], ncol(x), nsamp)
    }

    valmat <- unlist(lapply(z, function(x) x$valmat))
    dim(valmat) <- c(length(scheme@valname), nk, nsamp)
    dimnames(valmat)[[1]] <- scheme@valname
    dimnames(valmat)[[2]] <- k

    if(verbose) cat("\n")

    new("resampleFlexclust", k=as.integer(k), centers1=cent1, centers2=cent2,
        cluster1=clust1, cluster2=clust2, index1=index1, index2=index2,
        indextest=indextest, validation=valmat, call=MYCALL)
}


###**********************************************************

setMethod("show", "resampleFlexclust",
function(object){
    cat("An object of class", sQuote(class(object)),"\n\n")
    cat("Call:\n")
    print(object@call)
    cat("\nNumber of sample pairs:", length(object@index1),"\n")
})

setMethod("summary", "resampleFlexclust",
function(object){
    cat("Call:\n")
    print(object@call)
    cat("\n")
    ## only one number of clusters
    if(dim(object@validation)[2]==1){
        summary(t(object@validation[,1,]))
    } else {
        for(n in dimnames(object@validation)[[1]]){
            cat("Summary of", n, "\n")
            print(summary(t(object@validation[n,,])))
        }
    }
})

setMethod("plot", signature("resampleFlexclust","missing"),
function(x, y, ...){
    boxplot(x, ...)
})

setMethod("boxplot", "resampleFlexclust", 
function(x, which=1, ylab=NULL, ...){
    ## only one number of clusters
    if(dim(x@validation)[2]==1){
        boxplot(t(x@validation[,1,]), ylab=ylab, ...)
    } else {
        if(is.null(ylab)){
            if(is.character(which)) ylab <- which
            else ylab <- dimnames(x@validation)[[1]][which]
        }
        boxplot(as.data.frame(t(x@validation[which,,])), ylab=ylab, ...)
    }
})

setMethod("densityplot", "resampleFlexclust",
function(x, data, which=1, ...){
    ## only one number of clusters
    if(dim(x@validation)[2]==1){
        val <- t(x@validation[,1,])
    } else {
        val <- t(x@validation[which,,])
    }
    VAL <- as.vector(val)
    k <- rep(colnames(val), rep(nrow(val), ncol(val)))
    k <- factor(k, levels=colnames(val))

    densityplot(~VAL|k, as.table=TRUE, to=1, ...)
})


###**********************************************************

bsample <- function(object, ...)
{
    data <- getData(object)
    n <- nrow(data)
    d <- ncol(data)

    r <- matrix(rnorm(n*d), nrow=n)
    rs <- sqrt(rowSums(r^2))
    r <- r/rs

    e <- object@cldist[,1]
    e <- sample(e)
    r <- r*e

    ## transpose to fill in column-wise
    r <- t(r)
    for(kk in 1:object@k)
        r[,object@cluster==kk] <- r[,object@cluster==kk]+object@centers[kk,]

    t(r)
}
    
distsum <- function(object, data)
{
    distmat <- object@family@dist(data, object@centers)
    sum(distmat[cbind(1:nrow(data), object@cluster)])
}

setGeneric("paraBoot",
function(object, ...) standardGeneric("paraBoot"))

setMethod("paraBoot", signature(object="kcca"),
function(object, object2, nsamp=1000)
{
    x <- getData(object)
    stopifnot(all.equal(x, getData(object)))
    stopifnot(object@k < object2@k)

    ## abs should not be necessary, just to be on the safe side
    t0 <- abs(distsum(object, x) - distsum(object2, x))
    t1 <- double(nsamp)
    
    for(n in 1:nsamp)
    {
        bs1 <- bsample(object)
        t1[n] <-  abs(distsum(object, bs1) - distsum(object2, bs1))
    }
    mean(t1>=t0)
})
        
setMethod("paraBoot", signature(object="stepFlexclust"),
function(object, method=c("linear", "pairwise"), padjust=NULL, ...)
{
    method <- match.arg(method)
    k <- sort(object@k)
    nk <- length(k)
    if(method=="linear")
    {
        z <- double(nk-1)
        for(m in 1:(nk-1)){
            z[m] <- paraBoot(object[[m]], object[[m+1]], ...)
        }
        names(z) <- paste(k[-nk],k[-1],sep=":")
    }
    else {
        z <- matrix(NA, nrow=nk-1, ncol=nk-1)
        for(m in 1:(nk-1)){
            for(n in (m+1):nk){
                z[m,n-1] <- paraBoot(object[[m]], object[[n]], ...)
            }
        }
        colnames(z) <- k[-1]
        rownames(z) <- k[-nk]
    }
    if(!is.null(padjust))
        z[!is.na(z)] <- p.adjust(z[!is.na(z)], method=padjust)
    z
})
        
###**********************************************************

slswFlexclust <- function(x, object, ...)
{
    scheme <- clusterwiseScheme(object)
    resampleFlexclust(x=x, k=object@k, scheme=scheme, ...)
}
