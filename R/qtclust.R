qtclust <- function(x, radius, family=kccaFamily("kmeans"),
                    control=NULL)
{
    MYCALL <- match.call()
    control <- as(control, "flexclustControl")

    x <- as(x, "matrix")
    x <- family@preproc(x)

    cluster <- rep(NA, nrow(x))
    IND=1:nrow(x)
    ntry <- control@ntry

    k <- 1
    iter <- 0
    while(any( ok <- is.na(cluster) )){
        iter <- iter+1

        if(sum(ok)==1){
            cluster[ok] <- 0
            break
        }
        
        ntry <- min(ntry, sum(ok))

        ## first try to find a candidate centroid
        nfound <- 0
        for(n in 1:ntry){
            if(ntry==sum(ok))
               m = which(ok)[n]
            else
                m = sample(which(ok), 1)
            
            d = family@dist(x[ok,,drop=FALSE], x[m,,drop=FALSE])

            if(sum(d<=radius)>nfound){
                dok <- d<=radius
                nfound <- sum(dok)
            }
        }

        ## now try to add as many points as possible to the cluster
        ## without increasing the diameter.
        if(sum(dok)>1){
            olddok <- rep(FALSE, length(dok))
            while(any(olddok != dok)){
                olddok <- dok
                ## x[ok,][dok,] gives current cluster members
                d <- family@dist(x[ok,,drop=FALSE],
                                 x[ok,,drop=FALSE][dok,,drop=FALSE])

                ## create a set of candidates for addition, then
                ## iteratively remove all violating the diameter for
                ## the larger set
                dok <- d[cbind(1:nrow(d), max.col(d))] <= 2*radius
                cand <- which(dok != olddok)
                if(length(cand)>1){
                    d <- family@dist(x[ok,,drop=FALSE][cand,,drop=FALSE],
                                     x[ok,,drop=FALSE][cand,,drop=FALSE])
                    for(n in 1:length(cand)){
                        if(max(d[n,1:(n-1)])> 2*radius){
                            dok[cand[n]] <- FALSE
                            ## don't consider this point in the
                            ## following iterations
                            d[,n] <- 0
                        }
                    }
                }
            }
        }

        if(sum(dok)>control@min.size){
            cluster[ok][dok] <- k
            k <- k+1
        }
        else{
           cluster[ok][dok] <- 0
       }
    }

    ok <- cluster>0
    if(!any(ok)){
        stop("Could not find a valid clustering, try again with different radius")
    }

    cluster[!ok] <- NA
    cluster <- order(table(cluster), decreasing=TRUE)[cluster]
    centers <- family@allcent(x[ok,,drop=FALSE], cluster[ok])

    z <- newKccaObject(x, family, centers, converged=TRUE, call=MYCALL,
                       iter=as.integer(iter))

    ## newKccaObject assigns to closest centroid, we do not want this
    ## here -> reassign to original clusters.
    ## <FIXME>
    ##   the plot method will still use the KCCA assignment, not
    ##   sure if this is bug or feature
    ## </FIXME>
    z@cluster <- cluster
    z
}
