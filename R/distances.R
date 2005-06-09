distEuclidean <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        z[,k] <- sqrt( colSums((t(x) - centers[k,])^2) )
    }
    z
}

distManhattan <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        z[,k] <- colSums(abs(t(x) - centers[k,]))
    }
    z
}

distMax <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        z[,k] <- apply(abs(t(x) - centers[k,]), 2, max)
    }
    z
}

distJaccard <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    xc <- x %*% t(centers)
    nenner <-
        matrix(rowSums(x), nrow=nrow(x), ncol=nrow(centers)) +
            matrix(rowSums(centers), nrow=nrow(x), ncol=nrow(centers),
                   byrow=TRUE) - xc

    1 - xc/nenner
}

centOptim <- function(x, dist)
{
    foo <- function(p)
        sum(dist(x, matrix(p, nrow=1)))

    optim(colMeans(x), foo)$par
}

centOptim01 <- function(x, dist)
{
    foo <- function(p)
        sum(dist(x, matrix(p, nrow=1)))

    optim(colMeans(x), foo, lower=0, upper=1, method="L-BFGS-B")$par
}

###**********************************************************

distAngle <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        if(any(is.na(centers[k,])))
            z[,k] <- - Inf
        else
            z[,k] <- apply(x, 1, function(a) sum(a*centers[k,]))
    }
    1-z
}


centAngle <- function(x)
{
    z <- colMeans(x)
    z/sqrt(sum(z^2))
}

###**********************************************************

distCor <- function(x, centers)
{
   z <- matrix(0,nrow(x), ncol=nrow(centers))
   for(k in 1:nrow(centers)){
      z[,k] <- 1 - cor(t(x),centers[k,])
   }
   z
}   

