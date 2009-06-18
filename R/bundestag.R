##
## Copyright (C) 2008 Friedrich Leisch
## $Id: bundestag.R 4366 2009-06-18 13:43:44Z leisch $
##
bundestag <- function(year, second=TRUE, percent=TRUE, nazero=TRUE,
                      state=FALSE)
{
    year <- match.arg(as.character(year), c("2002", "2005"))
    tempenv <- new.env()
    if(year=="2002"){
        load(system.file("data/btw2002.RData", package="flexclust"),
             envir=tempenv)
        x <- get("btw2002", envir=tempenv)
    }
    else{
        load(system.file("data/btw2005.RData", package="flexclust"),
             envir=tempenv)
        x <- get("btw2005", envir=tempenv)
    }

    if(is.logical(state)){
        if(state) return(x$state)
    }
    else{
        y <- rep("other", nrow(x))
        ok <- grep(state, x$state)
        y[ok] <- as.character(x$state)[ok]
        return(as.factor(y))
    }

    if(second)
        p <- "2$"
    else
        p <- "1$"

    y <- x[,grep(p, colnames(x))]
    colnames(y) <- gsub(p, "", colnames(y))
    y <- as.matrix(y)
                        
    if(percent)
        y <- y/y[,"valid"]

    y <- y[,-grep("valid", colnames(y))]

    if(nazero) y[is.na(y)] <- 0

    y
}
