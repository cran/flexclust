%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: Nclus.Rd 4123 2008-09-25 11:27:34Z leisch $
%
\name{Nclus}
\alias{Nclus}
\title{Artificial Example with 4 Gaussians}
\description{
  A simple artificial regression example with 4 clusters, all of
  them having a Gaussian distribution.
}
\usage{
data(Nclus)
}
\details{
  The \code{Nclus} data set can be re-created by loading package
  \pkg{flexmix} and running \code{ExNclus(100)} 
  using \code{set.seed(2602)}. It has been saved as a data set for
  simplicity of examples only.
}
\examples{
data(Nclus)
cl <- cclust(Nclus, k=4, simple=FALSE, save.data=TRUE)
plot(cl)
}
\keyword{datasets}
