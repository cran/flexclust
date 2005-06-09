\name{flexclustControl-class}
\docType{class}
\alias{flexclustControl-class}
\alias{coerce,list,flexclustControl-method}
\alias{coerce,NULL,flexclustControl-method}
\alias{cclustControl-class}
\alias{coerce,list,cclustControl-method}
\alias{coerce,NULL,cclustControl-method}
\title{Classes "flexclustControl" and "cclustControl"}
\description{Hyperparameters for cluster algorithms.}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("flexclustControl", ...)}. In addition, named lists can be
  coerced to \code{flexclustControl}
  objects, names are completed if unique (see examples).
}

\section{Slots}{
  Objects of class \code{flexclustControl} have the following slots:
  \describe{
    \item{\code{iter.max}:}{Maximum number of iterations.}

    \item{\code{tolerance}:}{The algorithm is stopped when the
      (relative) change of the optimization criterium is smaller than
      \code{tolerance}.}
    
    \item{\code{verbose}:}{If a positive integer, then progress is
      reported every \code{verbose} iterations. If 0,
      no output is generated during model fitting.}

    \item{\code{classify}:}{Character string, one of \code{"auto"},
      \code{"weighted"}, \code{"hard"} or \code{"simann"}.}
    
    \item{\code{gamma}:}{Gamma value for weighted hard competitve
      learning.}
    \item{\code{simann}:}{Parameters for simulated annealing
      optimization (only used when \code{classify="simann"}).}
    
    \item{\code{ntry}:}{Number of trials per iteration for QT clustering.}
    \item{\code{min.size}:}{Clusters smaller than this value are treated
      as outliers.}
  }
  Objects of class \code{cclustControl} inherit from
  \code{flexclustControl} and have the following additional slots:
  \describe{
    \item{\code{method}:}{Learning rate for hard competitive learning,
      one of \code{"polynomial"} or \code{"exponential"}.}
    \item{\code{pol.rate}:}{Positive number for polynomial learning rate
      of form \eqn{1/iter^par}.}
    \item{\code{exp.rate}}{Vector of length 2 with paramters for
      exponential learning rate of form
      \eqn{par1*(par2/par1)^(iter/iter.max)}}.
    \item{\code{ng.rate}:}{Vector of length 3 with parameters for neural
      gas.}
  }
}
\author{Friedrich Leisch}
\seealso{
  \code{\link{kcca}}, \code{\link{cclust}}
}
\keyword{classes}
\examples{
## have a look at the defaults
new("flexclustControl")

## corce a list
mycont = list(iter=500, tol=0.001, class="w")
as(mycont, "flexclustControl")

## some additional slots
as(mycont, "cclustControl")
}
