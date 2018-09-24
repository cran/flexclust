#
#  Copyright (C) 2012 Friedrich Leisch
#  $Id: relabel.R 248 2018-04-30 09:36:44Z gruen $
#

dorelabel <- function(object, perm)
{
    po <- order(perm)

    object@cluster <- po[object@cluster]
    object@centers <- object@centers[perm,]
    object@clusinfo <- object@clusinfo[perm,,drop=FALSE]
    rownames(object@clusinfo) <- 1:object@k

    if(is(object, "kcca")){
        object@second <- po[object@second]
        object@clsim <- object@clsim[perm,perm]
    }
    
    object
}

setMethod("relabel", signature(object="kccasimple", by="character"),
function(object, by, which=NULL, ...)
{
    by <- match.arg(by, c("mean", "median", "variable", "manual",
                          "centers", "shadow", "symmshadow"))
    MYCALL <- match.call()

    # browser()
    if(by=="mean"){
        perm <- order(rowSums(object@centers))
    }
    else if(by=="median"){
        perm <- order(apply(object@centers, 1, median))
    }
    else if(by=="manual"){
        perm <- which
    }
    else if(by=="variable"){
        which <- which[1]
        perm <- order(object@centers[,which])
    }
    else{
        if(by == "centers"){
            d <- as.dist(1-clusterSim(object, method="centers"))
        }
        else if(by == "symmshadow"){
            d <- as.dist(1-clusterSim(object, method="shadow", symmetric=TRUE))
        }
        if(by == "shadow"){            
            d <- 1 - clusterSim(object, method="shadow", symmetric=FALSE)
        }
        perm <- as.integer(seriation::seriate(d, ...)[[1]])
    }
    object <- dorelabel(object, perm=perm)
    object@call <- MYCALL
    object
})


setMethod("relabel", signature(object="kccasimple", by="missing"),
function(object, by, ...)
{   
    MYCALL <- match.call()
    object <- relabel(object, by="mean", ...)
    object@call <- MYCALL
    object
})

setMethod("relabel", signature(object="kccasimple", by="integer"),
function(object, by, ...)
{
    MYCALL <- match.call()
    if(!all(sort(by) == 1:object@k))
        stop("if integer, ", sQuote("by"),
             " must be a permuatation of the numbers 1 to ", object@k)

    object <- dorelabel(object, by)
    object@call <- MYCALL
    object
})

###**********************************************************

setMethod("relabel", signature(object="stepFlexclust", by="character"),
function(object, by="series", ...)
{
    MYCALL <- match.call()

    if(by=="series")
    {
        object@models[[1]] <- relabel(object@models[[1]])
        for(k in 2:length(object@models))
        {
            TAB <- table(object@models[[k-1]]@cluster,
                         object@models[[k]]@cluster)
            ser <- mincross(TAB)
            object@models[[k]] <- dorelabel(object@models[[k]], ser)
        }
    }
    else
    {
        for(k in seq_along(object@models))
            object@models[[k]] <- relabel(object@models[[k]], by=by, ...)
    }
    object@call <- MYCALL
    object
})

mincross <- function(x)
{
    s <- clue::solve_LSAP(x, maximum=TRUE)

    ## columns not assigned to a cluster in the smaller solution (=rows)
    mm <- which(!(1:ncol(x)) %in% s)
    ## row with maximum in respective columns
    xm <- apply(x[,mm,drop=FALSE], 2, max)

    ## note: if dim(x) (k-1) times k then length(mm)=1
    ## while only needed when there are gaps in the number of clusters
    while(length(mm))
    {
        wm <- which.max(xm)
        m <- mm[wm]
        xm <- xm[-wm]
        mm <- mm[-wm]

        nr <- length(s)
        wohin <- which.max(x[,m])
        ## if whohin<nr insert right of target
        if(wohin<nr){
            s <- c(s[1:wohin], m, s[(wohin+1):nr])
        }
        else {  #else left
            s <- c(s[1:(nr-1)], m, s[nr])
        }
    }
    s
}

setMethod("relabel", signature(object="stepFlexclust", by="missing"),
function(object, by, ...)
{   
    MYCALL <- match.call()
    object <- relabel(object, by="series", ...)
    object@call <- MYCALL
    object
})

###**********************************************************

transitions <- function(object, cex=1.5, lwd=50, min=0.01,
                        enttext=FALSE, entcol=TRUE, slsa=FALSE,
                        power=NULL, hue=NULL)
{
    if(is.null(power)) power <- ifelse(slsa, 2, 1)
    if(is.null(hue)) hue <- ifelse(slsa, 0.4, 1)
    
    k <- object@k
    plot(0, xlim=c(1,length(k)), ylim=c(1,max(k)),
         xlab="", ylab="", type="n", axes=FALSE)

    rad <- (ceiling(log10(max(k))) + 2)*cex

    haveent <- is.list(entcol)
    if(haveent){
        MENT <- range(unlist(entcol))
        haveent <- TRUE
        ent2 <- entcol
        ent2 <- lapply(ent2, function(x) (x-MENT[1])/diff(MENT))
        entcol <- TRUE
        hue <- 0.65
        hue2 <- 1
    }

    ent <- list()

    bg <- hsv(1,0,0.85)
    for(n in 1:(length(k)-1)){
        x1 <- rep(n, k[n])
        y1 <- rev(seq(1, max(k), length=k[n]))
        x2 <- rep(n+1, k[n+1])
        y2 <- rev(seq(1, max(k), length=k[n+1]))

        z <- table(object@models[[n]]@cluster,
                   object@models[[n+1]]@cluster)

        z <- z/sum(z)
        
        ze <- t(t(z)/colSums(z))
        ## log(k) is entropy for discrete uniform with k levels
        ent[[n]] <- apply(ze, 2, entropy)/log(nrow(z))
        bgl <- rep("black",length(x2))

        if(slsa){
            if(entcol){
                if(haveent){
                    bg <- hsv(h=hue2, s=ent2[[n]]^power, v=.85)
                } else {
                    if(n>1)
                        bg <- hsv(h=hue, s=(1-ent[[n-1]])^power, v=.85)
                }
                bgl <- hsv(h=hue, s=(1-ent[[n]])^power, v=.85)
            }
        }
        else{
            if(enttext) text(x1+0.75, y2, round(ent[[n]], 3),
                             col="red", cex=cex)
            if(haveent){
                bg <- hsv(h=hue2, s=ent2[[n]]^power, v=.85)
            } else {
                if(n>1 & entcol){
                    bg <- hsv(h=hue, s=(ent[[n-1]])^power, v=.85)
                }
            }
        }
        for(i in seq_along(x1)){
            for(j in seq_along(x2)){
                if(z[i,j]>min)
                    segments(x1[i],y1[i],x2[j],y2[j],
                             lwd=lwd*z[i,j], col=bgl[j])
            }
        }

        points(x1, y1, pch=21, cex=rad, bg=bg)
        text(x1, y1, labels=1:k[n], cex=cex)
    }
    n <- n+1
    x1 <- rep(n, k[n])
    y1 <- rev(seq(1, max(k), length=k[n]))
    if(slsa){
        if(haveent){
            bg <- hsv(h=hue2, s=ent2[[n]]^power, v=.85)
        } else {
            if(n>1 & entcol)
                bg <- hsv(h=hue, s=(1-ent[[n-1]])^power, v=.85)
        }
    }
    else{
        if(enttext) text(x1+0.75, y2, round(ent[[n]], 3),
                         col="red", cex=cex)
        if(haveent){
            bg <- hsv(h=hue2, s=ent2[[n]]^power, v=.85)
        } else {
            if(n>1 & entcol)
                bg <- hsv(h=hue, s=(ent[[n-1]])^power, v=.85)
        }
    }
    points(x1, y1, pch=21, cex=rad, bg=bg)
    text(x1, y1, labels=1:k[n], cex=cex)
    invisible(ent)
}

slsaplot <- function(object, nodecol=NULL, ...)
{
    entcol=TRUE
    if(!is.null(nodecol))
        entcol <- clusterValues(object, nodecol)
    transitions(object, ..., entcol=entcol, slsa=TRUE)
}

entropy <- function(p)
{
    p <- p[p>0]
    -sum(p*log(p))
}

clusterValues <- function(object, x){
    z <- list()
    k <- 1
    for(cl in object@models){
        z[[k]] <- tapply(x, clusters(cl), mean, na.rm=TRUE)
        k <- k+1
    }
    names(z) <- names(object@models)
    z
}
