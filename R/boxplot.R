#
#  Copyright (C) 2013 Friedrich Leisch
#  $Id: boxplot.R 222 2017-03-03 16:29:43Z leisch $
#

setMethod("bwplot", "kccasimple",
function(x, data, xlab="", strip.labels=NULL, strip.prefix="Cluster ",
         col=NULL, shade=!is.null(shadefun), shadefun=NULL, byvar=FALSE, ...)
{
    if(is.character(shadefun))
        shadefun <- get(shadefun, mode="function")
    
    if(shade && is.null(shadefun)){
        shadefun <- medianInside
    }
    
    if(is.null(strip.labels)){
        strip.labels <- paste(strip.prefix, 1:x@k, sep="")
        if(!byvar){
            SIZE <- info(x, "size")
            strip.labels <-
            paste(strip.labels, ": ", SIZE, " (",
                  round(100 * SIZE/sum(SIZE)), "%)", sep="")
        }
    }

    if(missing(data)) data <- NULL

    if(is.null(col))
        col <- flxColors(color="full")

    col <- rep(col, length=x@k)

    DF <- kcca2df(x, data)
    levels(DF$group) <- strip.labels

    panel <- createBWplotPanel(data=DF, col=col,
                               shade=shade, shadefun=shadefun,
                               byvar=byvar)

    if(byvar){
        bwplot(group~value|variable, data=DF,
               as.table=TRUE, xlab=xlab, panel=panel, ...)
    } else {
        bwplot(variable~value|group, data=DF,
               as.table=TRUE, xlab=xlab, panel=panel, ...)
    }
})

###**********************************************************

createBWplotPanel <- function(data, col, shade, shadefun, byvar, ...)
{
    KKK <- 1
    KKKplus <- function() KKK <<- KKK+1
            
    grey <- flxColors(color="dark", grey=TRUE)
    nvar <- ifelse(byvar, nlevels(data$group),
                   nlevels(data$variable))

    brect1 <- brect2 <- trellis.par.get("box.rectangle")
    bumbr1 <- bumbr2 <- trellis.par.get("box.umbrella")
    
    bumbr2$col <- brect2$col <- grey
    bumbr1$lty <- bumbr2$lty <- 1

    if(shade){
        OK <- shadefun(data, ...)
    } else {
        OK <- matrix(TRUE, nrow=length(levels(data$variable)),
                     ncol=length(levels(data$group)))
    }

    mypanel <- function(x, y, ...)
    {
        if(shade){
            COL <- rep("white", nvar)
            if(byvar)
                COL[OK[KKK,]] <- col[OK[KKK,]]
            else
                COL[OK[,KKK]] <- col[KKK]
        }
        else{
            COL <- if(byvar) col else col[KKK]
        }

        trellis.par.set("box.rectangle", brect2)
        trellis.par.set("box.umbrella", bumbr2)
        if(!byvar)
            panel.bwplot(data$value, data$variable, box.ratio=0.3, col=grey)
        
        trellis.par.set("box.rectangle", brect1)
        trellis.par.set("box.umbrella", bumbr1)
        panel.bwplot(x, y, fill=COL, ...)

        KKKplus()
    }
    return(mypanel)
}

medianInside <- function(data, ...)
{
    ## median of population
    mp <- as.vector(tapply(data$value, list(data$variable), median))

    ## quartiles of groups
    qcl <- tapply(data$value,
                  list(data$variable, data$group), quantile, probs=0.25)
    qcu <- tapply(data$value, list(data$variable, data$group),
                  quantile, probs=0.75)

    (qcl>mp) | (qcu<mp)
}

boxOverlap <- function(data, ...)
{
    ## quartiles of population
    qpl <- as.vector(tapply(data$value, list(data$variable),
                            quantile, probs=0.25))
    qpu <- as.vector(tapply(data$value, list(data$variable),
                            quantile, probs=0.75))

    ## quartiles of groups
    qcl <- tapply(data$value, list(data$variable, data$group),
                  quantile, probs=0.25)
    qcu <- tapply(data$value, list(data$variable, data$group),
                  quantile, probs=0.75)

    (qcl>qpu) | (qcu<qpl)
}

###**********************************************************

kcca2df <- function(object, data=NULL)
{
    if(is.null(data)){
        data <- getData(object, error=TRUE)
        cluster <- clusters(object)
    }
    else
    {
        data <- as.matrix(data)
        cluster <- predict(object, data)
    }
    if(is.null(colnames(data)))
        colnames(data) <- paste("V", 1:object@k, sep="")

    matrix2df(data, cluster)
}

matrix2df <- function(data, group)
{
    data <- as.matrix(data)
    group <- as.factor(group)
    
    DF <- data.frame(value=as.vector(data),
                     variable=factor(rep(colnames(data),
                                         rep(nrow(data), ncol(data))),
                                     levels=colnames(data)))
    cbind(DF, group=group)
}
    
###**********************************************************

groupBWplot <- function(x, g, alpha=0.05, correct="holm", xlab="", col=NULL,
                        shade=!is.null(shadefun), shadefun=NULL,
                        strip.prefix="", strip.labels=NULL, which=NULL,
                        byvar=FALSE, ...)
{
    ## body very similar to bwplot,kcca-method -> fix bugs in both
    
    x <- as.matrix(x)
    g <- as.factor(g)
    
    if(!is.null(which))
        x <- x[,which,drop=FALSE]

    if(is.character(shadefun))
        shadefun <- get(shadefun, mode="function")
    
    if(shade && is.null(shadefun)){
        shadefun <- kruskalTest
    }

    if(is.null(col))
        col <- flxColors(color="full")[1]
    
    col <- rep(col, length=nlevels(g))

    DF <- matrix2df(x, g)

    panel <- createBWplotPanel(data=DF, col=col,
                               shade=shade, shadefun=shadefun,
                               alpha=alpha, correct=correct, byvar=byvar)

    if(byvar){
        bwplot(group~value|variable, data=DF,
               as.table=TRUE, xlab=xlab, panel=panel, ...)
    } else {
        bwplot(variable~value|group, data=DF,
               as.table=TRUE, xlab=xlab, panel=panel, ...)
    }
}
    
kruskalTest <- function(data, alpha, correct)
{
    nv <- nlevels(data$variable)
    pval <- double(nv)
    res <- matrix(FALSE, nrow=nv, ncol=nlevels(data$group))
    for(n in 1:nv)
    {
        pval[n] <- kruskal.test(value~group, data=data,
                                subset=data$variable==levels(data$variable)[n])$p.value
    }
    pval <- p.adjust(pval, method=correct)
    res[pval<alpha,] <- TRUE
    res
}
