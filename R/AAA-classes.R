#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: AAA-classes.R 3017 2006-10-02 12:45:13Z leisch $
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
                        call="call",
                        control="flexclustControl",
                        data="ModelEnv"))

setClass("kccasimple",
         contains="flexclust",
         representation(centers="ANY",
                        family="kccaFamily",
                        cldist="matrix"))       

setClass("kcca",
         contains="kccasimple",
         representation(second="integer",
                        xrange="ANY",           # range of all data
                        xcent="ANY",            # centroid of all data
                        totaldist="numeric",    # total dist data<->xcent
                        clsim="matrix"))

                        
                        
###**********************************************************

setClass("stepFlexclust",
         representation(models="list",
                        k="integer",
                        nrep="integer",
                        call="call",
                        xcent="ANY",            # centroid of all data
                        totaldist="numeric"     # total dist data<->xcent
                        ))


         
