setMethod("plot", signature(x="kcca",y="missing"),
function(x, y, which=1:2, project=NULL,
         data=NULL, points=!is.null(data), hull=!is.null(data),
         number = TRUE, simlines=TRUE,
         lwd=1, maxlwd=8*lwd, cex=1.5, numcol="black", nodes=16,
         add=FALSE, xlab="", ylab="", xlim = NULL,
         ylim = NULL, pch=NULL, col=NULL, ...)
{
    if(length(which)!=2)
        stop(sQuote("which"), " must have length 2")

    if(is.null(project))
        project <- function(x) x
    else{
        if(!is(project, "function")){
            probject <- project
            project <- function(x) predict(probject, x)
        }
    }    
    centers <- project(x@centers)[,which]
    xrange <- apply(centers, 2, range)

    if(is(hull, "function")){
        HULLFUN <- hull
        hull <- !is.null(data)
    }
    else if(is(hull, "character")){
        HULLFUN <- get(hull, mode="function")
        hull <- !is.null(data)
    }        
    else{
        ## change to clusterHulls when published
        HULLFUN <- clusterEllipses
    }

    if(is.null(col)) col <- rep(1:8, length=x@k)

    if(!is.null(data)){
        distmat <- x@family@dist(data, x@centers)
        cluster <- x@family@cluster(distmat=distmat)
        cldist <- distmat[cbind(1:nrow(data), cluster)]

        data <- project(data)[,which]
        xrange <- apply(data, 2, range)
    }

    if(points){
        if(is.null(pch)){
            pch <- (cluster %% 10)
            pch[pch==0] <- 10
        }
        else if(length(pch)!=nrow(data)){
            pch <- rep(pch, length=x@k)
            pch <- pch[cluster]
        }
    }
    
    if(!add){
        if(points && !is.null(data))
            plot(data, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                 pch=pch, col=col[cluster], ...)
        else{
            if(is.null(xlim)) xlim <- sort(xrange[,1])
            if(is.null(ylim)) ylim <- sort(xrange[,2])
            plot(0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
        }
    }
    else{
        if(points)
            points(data, pch=pch, ...)
    }

    if(hull && !is.null(data)) HULLFUN(data, cluster, cldist, col)
    
    if(simlines){
        sim <- maxlwd*(x@clsim+t(x@clsim))/2

        for(k in 1:(x@k-1)){
            for(m in k:x@k)
                if(sim[k,m]>0)
                    lines(centers[c(k,m),], lwd=sim[k,m])
        }
    }

    if(number){
        rad <- ceiling(log10(x@k)) + 1.5
        points(centers, pch=21, cex=rad*cex, bg="white")
        text(centers, labels=1:x@k, cex=cex, col=numcol)
    }
    else
        points(centers, pch=nodes)

    invisible()
})

clusterHulls <- function(data, cluster, dist, col)
{
    col <- rep(col, length=max(cluster))

    for(k in 1:max(cluster)){
        ok <- cluster==k
        if(length(ok)>3){
            ok1 <- ok & (dist < median(dist[ok]))
            if(length(ok1)>3){
                hpts <- chull(data[ok1,])
                polygon(data[ok1,][c(hpts, hpts[1]),],
                        border=col[k], lwd=2*par("lwd"))
            }
            ok1 <- ok & (dist < 2.5*median(dist[ok]))
            if(length(ok1)>3){
                hpts <- chull(data[ok1,])
                polygon(data[ok1,][c(hpts, hpts[1]),], border=col[k], lty=2)
            }
        }
    }
}

clusterEllipses <- function(data, cluster, dist, col)
{
    require("ellipse")
    col <- rep(col, length=max(cluster))
    for(k in 1:max(cluster)){
        ok <- cluster==k
        if(length(ok)>3){
            lines(ellipse(cov(data[ok,]), centre=colMeans(data[ok,])),
                  col=col[k])
        }
    }
}
    

###**********************************************************

setMethod("image", signature(x="kcca"),
function(x, which = 1:2, npoints = 100,
         xlab = "", ylab = "", fastcol = TRUE, col=NULL, ...) 
{
    if(length(which)!=2)
        stop(sQuote(which), "must be a vector of length 2")
    
    X <- seq(x@xrange[1,which[1]],
             x@xrange[2,which[1]], length=npoints)
    Y <- seq(x@xrange[1,which[2]],
             x@xrange[2,which[2]], length=npoints)
    Y1 <- rep(Y, rep.int(length(X), length(Y)))
    X1 <- rep(X, times = ceiling(length(Y1)/length(X)))
    
    neib1 <- neighbours(x)

    if(fastcol){
        Z <- x@family@cluster(cbind(X1,Y1), x@centers[,which], n=1)
        neib2 <- neib1
    }
    else{
        Z <- x@family@cluster(cbind(X1,Y1), x@centers[,which], n=2)
        neib2 <- neighbours(Z)
        Z <- Z[[1]]
    }
    
    image(X, Y, matrix(Z, nrow=length(X)),
          col=neighbourColors(neib2, col=col), xlab=xlab, ylab=ylab)

    sim <- (x@clsim+t(x@clsim))/2
    plot(x, which=which, add=TRUE, ...)
    
    invisible(list(x=X,y=Y,z=Z))
})



neighbourColors <- function(object, col=NULL)
{
    ## library(colorspace)
    ## x = hcl(seq(0, 360*7/8, length = 8), c=30)
    ## dput(x[c(1,3,5,7,2,4,6,8)])
    if(is.null(col))
        col <- c("#F9C3CD", "#D0D4A8", "#9DDDD5", "#D1CCF5",
                 "#EDCAB2", "#AFDCB8", "#ACD7ED", "#EFC4E8")
    
    icol <- seq(along=col)
    z <- rep(NA, length(object))
    names(z) <- 1:length(z)

    nmis <- 0
    while(any(is.na(z))){
        
        node <- min(which(is.na(z)))
        nodes <- c(node,object[[node]])
        
        for(n in nodes){
            if(is.na(z[n])){
                used <- na.omit(z[object[[n]]])
                if(length(used)>=length(col)){
                    nmis <- nmis+1
                    z[n] <- sample(icol, 1)
                }
                else if(length(used>0))
                    z[n] <- icol[-used][1]
                else
                    z[n] <- icol[1]
            }
        }
    }
    if(nmis>0)
        warning("colors were missing ", nmis, " times\n")
    
    col[z]
}

###**********************************************************

setMethod("barplot", "kcca",
function (height, bycluster = TRUE, oneplot = TRUE,
          data = NULL, FUN=colMeans, 
          main = deparse(substitute(height)), 
          which = 1:height@k,
          names.arg = NULL, xcent=TRUE, oma=par("oma"),
          ...) 
{
    object <- height
    opar <- par(c("mfrow", "oma", "mgp", "xpd", "oma"))
    on.exit(par(opar))
    par(xpd = NA)
    n <- length(which)

    if(is.null(data)){
        cluster <- object@cluster
        centers <- object@centers
        size <- object@size
        datacent <- object@xcent
    }
    else{
        cluster <- predict(object, data)
        centers <- matrix(NA, nrow=object@k, ncol=ncol(data))
        colnames(centers) <- colnames(data)
        for(k in 1:object@k){
            ok <- cluster==k
            if(any(ok))
                centers[k,] <- FUN(data[ok,])
        }
        size <- tabulate(cluster, nbins=object@k)
        datacent <- FUN(data)
    }
    
    ylim <- range(centers)
    
    if (is.null(names.arg)) 
        names.arg <- colnames(centers)
    if (bycluster) {
        if (oneplot) {
            if (n <= 3) {
                par(mfrow = c(n, 1), oma=oma)
            }
            else {
                par(mfrow = c(ceiling(n/2), 2), oma=oma)
            }
        }
        for (k in which) {
            mid <- barplot(centers[k, ], col = "grey", 
                names.arg = NULL, ylim = ylim, ...)
            if (!is.null(names.arg)) {
                text(mid + 0.3 * (mid[2] - mid[1]), 0, adj = 1, 
                  srt = 45, paste(names.arg, "  "))
            }
            if (xcent && !is.null(datacent)) {
                points(mid, datacent, pch = 16, col = "red")
                points(mid, datacent, type = "h", col = "red")
            }
            title(main = paste("Cluster ", k, ": ", size[k], 
                " points (", round(100 * size[k]/sum(size), 
                  2), "%)", sep = ""))
        }
    }
    else {
        a <- ceiling(sqrt(ncol(centers)))
        if (oneplot) {
            par(mfrow = c(a, ceiling(ncol(centers)/a)))
        }
        for (k in 1:ncol(centers)) {
            barplot(centers[which, k, drop = FALSE], col = "grey", 
                    ylim = ylim, ...)
            title(main = names.arg[k])
            if (xcent && !is.null(datacent)) {
                abline(h = datacent[k], col = "red")
            }
        }
    }
})


###**********************************************************

setMethod("plot", signature(x="stepFlexclust", y="missing"),
function(x, y, xlab=NULL, ylab=NULL, type=c("barplot", "lines"),
         ...)
{
    type <- match.arg(type)
    X <- x@K
    Y <- sapply(x@models, function(z) sum(z@withindist))

    if(is.null(xlab))
        xlab <- "number of clusters"
    if(is.null(ylab))
         ylab <- "sum of within cluster distances"

    
    if (type == "barplot") 
        barplot(Y, names = X, xlab=xlab, ylab = ylab, ...)
    else {
        plot(X, Y, type = "b", axes = FALSE,
             xlab = xlab, ylab = ylab, ...)
        axis(2)
        axis(1, at = X, labels = X)
    }
    invisible()
})
