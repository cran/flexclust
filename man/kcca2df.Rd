\name{kcca2df}
\alias{kcca2df}
\title{Convert Cluster Result to Data Frame}
\description{
  Convert object of class \code{"kcca"} to a data frame in long format.
}
\usage{
kcca2df(object, data)
}
\arguments{
  \item{object}{Object of class \code{"kcca"}.}
  \item{data}{Optional data if not saved in \code{object}.}
}
\value{
  A \code{data.frame} with columns \code{value}, \code{variable} and
  \code{group}.
}
\examples{
c.iris <- cclust(iris[,-5], 3, save.data=TRUE)
df.c.iris <- kcca2df(c.iris)
summary(df.c.iris)
densityplot(~value|variable+group, data=df.c.iris)
}
\keyword{cluster}
