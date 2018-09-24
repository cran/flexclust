\name{auto}
\alias{auto}
\docType{data}
\title{
  Automobile Customer Survey Data
}
\description{
  A German manufacturer of premium cars asked customers approximately 3
  months after a car purchase which characteristics of the car were most
  important for the decision to buy the car. The survey was done in 1983
  and the data set contains all responses without missing values.
}
\usage{data(auto)}
\format{
  A data frame with 793 observations on the following 46 variables.
  \describe{
    \item{\code{model}}{A factor with levels \code{A}, \code{B},
      \code{C}, or \code{D}; model bought by the customer.}
    \item{\code{gear}}{A factor with levels \code{4 gears}, \code{5
	econo}, \code{5 sport}, or \code{automatic}.}
    \item{\code{leasing}}{A logical vector, was leasing used to finance
      the car?}
    \item{\code{usage}}{A factor with levels \code{private}, \code{both}, \code{business}.}
    \item{\code{previous_model}}{A factor describing which type of car
      was owned directly before the purchase.}
    \item{\code{other_consider}}{A factor with levels \code{same manuf},
      \code{other manuf}, \code{both}, or \code{none}.}
    \item{\code{test_drive}}{A logical vector, did you do a test drive?}
    \item{\code{info_adv}}{A logical vector, was advertising an
      important source of information?}
    \item{\code{info_exp}}{A logical vector, was experience an
      important source of information?}
    \item{\code{info_rec}}{A logical vector, were recommendations an
      important source of information?}
    \item{\code{ch_clarity}}{A logical vector.}
    \item{\code{ch_economy}}{A logical vector.}
    \item{\code{ch_driving_properties}}{A logical vector.}
    \item{\code{ch_service}}{A logical vector.}
    \item{\code{ch_interior}}{A logical vector.}
    \item{\code{ch_quality}}{A logical vector.}
    \item{\code{ch_technology}}{A logical vector.}
    \item{\code{ch_model_continuity}}{A logical vector.}
    \item{\code{ch_comfort}}{A logical vector.}
    \item{\code{ch_reliability}}{A logical vector.}
    \item{\code{ch_handling}}{A logical vector.}
    \item{\code{ch_reputation}}{A logical vector.}
    \item{\code{ch_concept}}{A logical vector.}
    \item{\code{ch_character}}{A logical vector.}
    \item{\code{ch_power}}{A logical vector.}
    \item{\code{ch_resale_value}}{A logical vector.}
    \item{\code{ch_styling}}{A logical vector.}
    \item{\code{ch_safety}}{A logical vector.}
    \item{\code{ch_sporty}}{A logical vector.}
    \item{\code{ch_consumption}}{A logical vector.}
    \item{\code{ch_space}}{A logical vector.}
    \item{\code{satisfaction}}{A numeric vector describing overall
      satisfaction (1=very good, 10=very bad).}
    \item{\code{good1}}{Conception, styling, dimensions.}
    \item{\code{good2}}{Auto body.}
    \item{\code{good3}}{Driving and coupled axles.}
    \item{\code{good4}}{Engine.}
    \item{\code{good5}}{Electronics.}
    \item{\code{good6}}{Financing and customer service.}
    \item{\code{good7}}{Other.}
    \item{\code{sporty}}{What do you think about the balance of
      sportiness and comfort? (\code{good}, \code{more sport}, \code{more comfort}).}
    \item{\code{drive_char}}{Driving characteristis (\code{gentle} < \code{speedy} < \code{powerfull} < \code{extreme}).}
    \item{\code{tempo}}{Which average speed do you prefer on German
      Autobahn in km/h? (\code{< 130} < \code{130-150} < \code{150-180} < \code{> 180})}
    \item{\code{consumption}}{An ordered factor with levels \code{low} < \code{ok} < \code{high} < \code{too high}.}
    \item{\code{gender}}{A factor with levels \code{male} and \code{female}}
    \item{\code{occupation}}{A factor with levels \code{self-employed},
      \code{freelance}, and \code{employee}.}
    \item{\code{household}}{Size of household, an ordered factor with levels \code{1-2} < \code{>=3}.}
  }
}
\source{
  The original German data are in the public domain and available
  from LMU Munich (doi:10.5282/ubm/data.14). The variable names and
  help page were translated to English and converted into Rd format
  by Friedrich Leisch.
}
\references{
  Open Data LMU (1983): Umfrage unter Kunden einer Automobilfirma,
  doi:10.5282/ubm/data.14
}
\examples{
data(auto)
summary(auto)
}
\keyword{datasets}
