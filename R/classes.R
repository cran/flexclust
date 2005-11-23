#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: classes.R 1779 2005-08-18 09:23:15Z leisch $
#

setClass("flexclustControl",
         representation(iter.max="numeric",
                        tolerance="numeric",
                        verbose="numeric",
                        classify="character",
                        gamma="numeric",          # for kcca
                        simann="numeric",         # for kcca
                        ntry="numeric",           # for qtclust
                        min.size="numeric"        # for qtclust
                        ),
         prototype(iter.max=200,
                   minprior=0.05,
                   tolerance=10e-7,
                   verbose=0,
                   classify="auto",
                   gamma=1,
                   simann=c(0.3, 0.95, 10),
                   ntry=5,
                   min.size=2))


setAs("list", "flexclustControl",
function(from, to){
    z <- list2object(from, to)
    z@classify <- match.arg(z@classify,
                            c("auto", "weighted", "hard", "simann"))
    z
})

setAs("NULL", "flexclustControl",
function(from, to){
    new(to)
})

###**********************************************************

setClass("cclustControl",
         contains="flexclustControl",
         representation(pol.rate="numeric",
                        exp.rate="numeric",
                        ng.rate="numeric",
                        method="character"),
         prototype(pol.rate=c(1,0),
                   exp.rate=c(0.1, 0.0001),
                   ng.rate=c(0.5, 0.005, 10, 0.01),
                   method="polynomial"))


setAs("list", "cclustControl",
function(from, to){
    z <- list2object(from, to)
    z@method <- match.arg(z@method,
                          c("polynomial", "exponential"))
    z
})

setAs("NULL", "cclustControl",
function(from, to){
    new(to)
})


###**********************************************************
###**********************************************************

setClass("kccaFamily",
         representation(name="character",
                        similarity="logical",
                        dist="function",
                        cent="function",
                        allcent="function",
                        wcent="function",
                        weighted="logical",
                        cluster="function",
                        preproc="function",
                        groupFun="function"),
         prototype(weighted=FALSE,
                   preproc=function(x) x))
                        


###**********************************************************

setClass("flexclust",
         representation(k="integer",
                        cluster="integer",
                        iter="integer",
                        converged="logical",
                        clusinfo="data.frame",
                        xrange="ANY",
                        call="call",
                        control="flexclustControl"))

setClass("kcca",
         contains="flexclust",
         representation(centers="ANY",
                        second="integer",
                        family="kccaFamily",
                        xcent="ANY",            # centroid of all data
                        totaldist="numeric",    # total dist data<->xcent
                        clsim="matrix",
                        cldist="matrix"))
                        
                        
###**********************************************************

setClass("stepFlexclust",
         representation(models="list",
                       K="integer",
                       nrep="integer",
                       call="call"))


         
