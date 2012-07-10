setGeneric("randIndex", function(x, y, correct=TRUE)
           standardGeneric("randIndex"))

setMethod("randIndex", signature(x="flexclust", y="flexclust"),
function(x, y, correct=TRUE){
    randIndex(table(clusters(x), clusters(y)), correct=correct)
})

setMethod("randIndex", signature(x="flexclust", y="integer"),
function(x, y, correct=TRUE){
    randIndex(table(clusters(x), y), correct=correct)
})

setMethod("randIndex", signature(x="integer", y="flexclust"),
function(x, y, correct=TRUE){
    randIndex(table(x, clusters(y)), correct=correct)
})

setMethod("randIndex", signature(x="integer", y="integer"),
function(x, y, correct=TRUE){
    randIndex(table(x, y), correct=correct)
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

