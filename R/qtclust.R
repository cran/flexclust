#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: qtclust.R 3887 2008-02-18 10:19:40Z leisch $
#

qtclust <- function(x, radius, family=kccaFamily("kmeans"),
                    control=NULL, simple=TRUE, save.data=FALSE)
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
    iter.width <- log10(nrow(x))+2
    
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
                dok <- (d <= radius)
                dok2 <- (d <= 2*radius)
                nfound <- sum(dok)
            }
        }
        
        ## We do not need to consider the following points
        ok[ok][!dok2] <- FALSE
        dok <- dok[dok2]

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
                dok[d[cbind(1:nrow(d), max.col(d))] <= 2*radius] <- TRUE
                cand <- which(dok != olddok)
                candok <- rep(FALSE, length(cand))
                candok[1] <- TRUE
                
                if(length(cand)>1){
                    for(n in 2:length(cand)){
                        d <- family@dist(x[ok,,drop=FALSE][cand[candok],,drop=FALSE],
                                         x[ok,,drop=FALSE][n,,drop=FALSE])

                        if(max(d) <= 2*radius){
                            candok[n] <- TRUE
                        }
                        else{
                            dok[cand[n]] <- FALSE
                        }
                    }
                }
            }
        }

        if(sum(dok)>=control@min.size){
            cluster[ok][dok] <- k
            k <- k+1
        }
        else{
            cluster[ok][dok] <- 0
        }

        if(control@verbose && (iter%%control@verbose==0))
            printIter(iter, sum(is.na(cluster)), "-- points remaining",
                      format="d", width=iter.width)

    }

    ok <- cluster>0
    if(!any(ok)){
        stop("Could not find a valid clustering, try again with different radius.")
    }
    if(all(cluster==1)){
        stop("All points in one cluster, try smaller radius.")
    }

    cluster[!ok] <- NA
    cluster <- order(table(cluster), decreasing=TRUE)[cluster]
    centers <- family@allcent(x[ok,,drop=FALSE], cluster[ok])

    if(simple)
    {
        z <- new("kccasimple",
                 k=nrow(centers),
                 cluster=cluster,
                 iter=as(iter, "integer"),
                 converged=TRUE,
                 call=MYCALL,
                 control=control,
                 centers=centers,
                 family=family,
                 clusinfo=data.frame(size=as.integer(table(cluster[ok]))))
    }
    else{
        z <- newKccaObject(x, family, centers, converged=TRUE, call=MYCALL,
                           iter=as.integer(iter))

        ## newKccaObject assigns to closest centroid, we do not want this
        ## here -> reassign to original clusters.
        ## <FIXME>
        ##   the plot method will still use the KCCA assignment, not
        ##   sure if this is bug or feature
        ## </FIXME>
        z@cluster <- cluster
    }

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=x)
    
    return(z)
}
