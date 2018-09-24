%
%  Copyright (C) 2012 Friedrich Leisch
%  $Id: priceFeature.Rd 252 2018-09-17 08:40:24Z gruen $
%
\name{priceFeature}
\alias{priceFeature}
\alias{plot.priceFeature}
\title{Artificial 2d Market Segment Data}
\description{
  Simple artificial 2-dimensional data to demonstrate clustering for market
  segmentation. One dimension is the hypothetical feature
  sophistication (or performance or quality, etc) of a product, the
  second dimension the price customers are willing to pay for the
  product.
}
\usage{
priceFeature(n, which=c("2clust", "3clust", "3clustold", "5clust",
                        "ellipse", "triangle", "circle", "square",
                        "largesmall"))
}
\arguments{
  \item{n}{Sample size.}
  \item{which}{Shape of data set.}
}
\references{
  Sara Dolnicar and Friedrich Leisch. Evaluation of structure and
  reproducibility of cluster solutions using the bootstrap. Marketing
  Letters, 21:83-101, 2010.
}
\keyword{datagen}
\examples{
plot(priceFeature(200, "2clust"))
plot(priceFeature(200, "3clust"))
plot(priceFeature(200, "3clustold"))
plot(priceFeature(200, "5clust"))
plot(priceFeature(200, "ell"))
plot(priceFeature(200, "tri"))
plot(priceFeature(200, "circ"))
plot(priceFeature(200, "square"))
plot(priceFeature(200, "largesmall"))
}
