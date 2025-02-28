#
#  Copyright (C) 2016 Friedrich Leisch
#
setMethod("histogram",signature(x="kccasimple", data="missing"),
function(x, data, xlab="", ...)
{
    data <- getData(x, error=TRUE)
    histogram(x, data=data, xlab=xlab, ...)
})

setMethod("histogram",signature(x="kccasimple", data="data.frame"),
function(x, data, xlab="", ...)
{
    histogram(x, data=data.matrix(data), xlab=xlab, ...)
})

setMethod("histogram",signature(x="kccasimple", data="matrix"),
function(x, data, xlab="Similarity", power=1, ...)
{
    dists <- x@family@dist(data, x@centers)
    dists <- dists[rowSums(is.na(dists)) != ncol(dists),]

    dists <- exp(-(dists)^power)
    sim <- sqrt(dists / rowSums(dists))

    ##sim <- 1 - exp(-dists)

    df <- data.frame(value=as.vector(sim))
    cluster <- rep(1:x@k, rep(nrow(data), x@k))
    df$cluster <- factor(paste("Cluster", cluster))
    
    histogram(~value|cluster, data=df, as.table=TRUE, xlab=xlab,
             ...)
})
