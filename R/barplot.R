#
#  Copyright (C) 2006 Friedrich Leisch
#  $Id: barplot.R 4259 2009-02-01 16:08:08Z leisch $
#

setMethod("barplot", "kccasimple",
function (height, bycluster = TRUE, oneplot = TRUE,
          data = NULL, FUN=colMeans, 
          main = deparse(substitute(height)), 
          which = 1:height@k,
          names.arg = NULL, oma=par("oma"),
          col=NULL, mcol="darkred", srt=45, horiz=FALSE, ...)
{
    object <- height
    opar <- par(c("mfrow", "oma", "mgp", "xpd", "oma"))
    on.exit(par(opar))
    n <- length(which)

    if(is.null(col))
        col <- LightColors
    
    col <- rep(col, length=object@k)

    if(is.null(data)){
        cluster <- object@cluster
        centers <- object@centers
        size <- info(object, "size")
        if(is(object, "kcca"))
            datacent <- object@xcent
        else
            datacent <- NULL
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
        par(xpd=NA)
        if (oneplot) {
            if (n <= 3) {
                par(mfrow = c(n, 1), oma=oma)
            }
            else {
                par(mfrow = c(ceiling(n/2), 2), oma=oma)
            }
        }
        for (k in which) {
            mid <- barplot(centers[k, ], col = col[k], 
                           names.arg = "", ylim = ylim, ...)
            if (!is.null(names.arg)){
                text(mid + 0.005 * min(srt, 90-srt) * (mid[2] - mid[1]),
                     ylim[1] - par("cxy")[2], adj = ifelse(srt==0, 0.5, 1),
                     srt = srt, paste(names.arg, "  "))
            }
            if (!is.null(mcol) && !is.null(datacent)) {
                points(mid, datacent, pch = 16, col = mcol)
                points(mid, datacent, type = "h", col = mcol)
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
            barplot(centers[which, k], col = col[which], 
                    ylim = ylim, xlab="", ...)
            title(main = names.arg[k])
            if (!is.null(mcol) && !is.null(datacent)) {
                abline(h = datacent[k], col = mcol)
            }
        }
    }
})

###**********************************************************

setMethod("barchart", "kccasimple",
function(x, data, xlab="", strip.prefix="Cluster ",
         col=NULL, mcol="darkred", panel=NULL, which=NULL, ...)
{
    if(is.null(col))
        col <- flxColors(color="light")

    if(is.null(which))
        which <- seq(1, ncol(x@centers))

    col <- rep(col, x@k)

    TAB <- x@centers[,which]
    ## sonst musz man die barplots von unten nach oben lesen
    TAB <- TAB[,ncol(TAB):1]
    xcent <- rev(x@xcent[which])


    SIZE <- info(x, "size")
    rownames(TAB) <-
        paste(strip.prefix, 1:nrow(TAB), ": ", SIZE, " (",
              round(100 * SIZE/sum(SIZE), 2), "%)", sep="")
    TAB <- as.data.frame(as.table(TAB))

    panel <- createBarchartPanel(m=xcent, col=col,
                                 mcol=mcol, which=which)
d
    barchart(Var2~Freq|Var1, data=TAB,
             panel=panel, as.table=TRUE,
             xlab=xlab, ...)
})

createBarchartPanel <- function(m, col, mcol)
{
    KKK <- 1
    KKKplus <- function() KKK <<- KKK+1

    ## sonst musz man die barplots von unten nach oben lesen

    mypanel <- function(x, y, shade=FALSE, diff=NULL, ...)
    {
        if(is.null(diff))
            diff <- c(max(m)/4, 0.5)
        else
            diff <- rep(diff, length=2)

        grey <- flxColors(color="dark", grey=TRUE)
        COL <- col[KKK]
        MCOL <- mcol
        BCOL <- "black"
            
        if(shade){
            COL <- rep("white", length(x))
            MCOL <- rep(grey, length=length(x))
            BCOL <- rep(grey, length=length(x))
            d1 <- abs(x-m) >= diff[1]
            d2 <- abs((x-m)/m) >= diff[2]

            COL[d1|d2] <- col[KKK]
            MCOL[d1|d2] <- mcol
            BCOL[d1|d2] <- "black"
        }
        
            
        grid.segments(x0=0, y0=1:length(x), x1=m, y1=1:length(x),
                      gp=gpar(col=MCOL),
                      default.units="native")

        panel.barchart(x, y, col=COL, border=BCOL, ...)
        grid.points(m, 1:length(x), pch=16,
                    size=unit(0.5, "char"), gp=gpar(col=MCOL))
        grid.segments(1, 1, 4, 4)
        KKKplus()
    }
    return(mypanel)
}
