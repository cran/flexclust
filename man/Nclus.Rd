%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: Nclus.Rd 1849 2005-10-10 06:15:57Z leisch $
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
plot(Nclus, col=rep(1:4, c(100, 100, 150, 200)))
}
\keyword{datasets}
