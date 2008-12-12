##
## Copyright (C) 2008 Friedrich Leisch
## $Id: bundestag.R 4207 2008-12-12 14:46:53Z leisch $
##
bundestag <- function(year, second=TRUE, percent=TRUE, nazero=TRUE, state=FALSE)
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

    if(state) return(x$state)

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
