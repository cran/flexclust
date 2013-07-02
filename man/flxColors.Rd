%
%  Copyright (C) 2005 Friedrich Leisch
%  $Id: flxColors.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{flxColors}
\alias{flxColors}
\title{Flexclust Color Palettes}
\description{
  Create and access palettes for the plot methods.
}
\usage{
  flxColors(n=1:8, color=c("full","medium", "light","dark"), grey=FALSE)
}
\arguments{
  \item{n}{Index number of color to return (1 to 8)}
  \item{color}{Type of color, see details.}
  \item{grey}{Return grey value corresponding to palette.}
}
\details{
  This function creates color palettes in HCL space for up to 8
  colors. All palettes have constant chroma and luminance, only the hue
  of the colors change within a palette.
  
  Palettes \code{"full"} and \code{"dark"} have the same luminance, and
  palettes \code{"medium"} and \code{"light"} have the same luminance.
}
\author{Friedrich Leisch}
\seealso{\code{\link{hcl}}}  
\keyword{color}
\examples{
opar <- par(c("mfrow", "mar", "xaxt"))
par(mfrow=c(2,2), mar=c(0,0,2,0), yaxt="n")

x <- rep(1, 8)

barplot(x, col = flxColors(color="full"), main="full")
barplot(x, col = flxColors(color="dark"), main="dark")
barplot(x, col = flxColors(color="medium"), main="medium")
barplot(x, col = flxColors(color="light"), main="light")

par(opar)
}
