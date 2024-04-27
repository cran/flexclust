\name{vacmot}
\alias{vacmot}
\alias{vacmot6}
\alias{vacmotdesc}
\docType{data}
\title{
  Vacation Motives of Australians
}
\description{
  
  In 2006 a sample of 1000 respondents representative for the adult Australian
  population was asked about their environmental
  behaviour when on vacation. In addition the survey also included a
  list of statements about vacation motives like "I want to rest and
  relax," "I use my holiday for the health and beauty of my body," and
  "Cultural offers and sights are a crucial factor.". Answers are binary
  ("applies", "does not apply").
}
\usage{data(vacmot)}
\format{
  Data frame \code{vacmot} has 1000 observations on 20 binary
  variables on travel motives. Data frame \code{vacmotdesc} has 1000
  observation on sociodemographic descriptor variables, mean moral
  obligation to protect the environment score, mean NEP score, and
  mean environmental behaviour score, see Dolnicar & Leisch
  (2008) for details.
  In addition integer vector \code{vacmot6} contains the 6
  cluster partition presented in Dolnicar & Leisch (2008).
}
\source{
The data set was collected by the Institute for Innovation in
Business and Social Research, University of Wollongong (NSW, Australia). 
}
\references{
  Sara Dolnicar and Friedrich Leisch. An investigation of tourists'
  patterns of obligation to protect the environment. Journal of Travel
  Research, 46:381-391, 2008.

  Sara Dolnicar and Friedrich Leisch. Using graphical statistics to
  better understand market segmentation solutions. International Journal
  of Market Research, 56(2):97-120, 2014.
}
\examples{
data(vacmot)
summary(vacmotdesc)
dotchart(sort(colMeans(vacmot)))

## reproduce Figure 6 from Dolnicar & Leisch (2008)
cl6 <- kcca(vacmot, k=vacmot6, control=list(iter=0))
barchart(cl6)
}
\keyword{datasets}
