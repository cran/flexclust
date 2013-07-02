%
%  Copyright (C) 2008 Friedrich Leisch
%  $Id: bundestag.Rd 3 2013-06-12 10:06:43Z leisch $
%
\name{bundestag}
\alias{bundestag}
\alias{btw2002}
\alias{btw2005}
\alias{btw2009}
\title{German Parliament Election Data}
\description{
  Results of the elections 2002, 2005 or 2009 for the German Bundestag, the
  first chamber of the German parliament.
}
\usage{
data(btw2002)
data(btw2005)
data(btw2009)
bundestag(year, second=TRUE, percent=TRUE, nazero=TRUE, state=FALSE)
}
\format{
  \code{btw200x} are data frames with 299 rows
  (corresponding to constituencies) and 17 columns. All columns except
  \code{state} are numeric.
  \describe{
    \item{state:}{factor, the 16 German federal states.}
    \item{eligible:}{number of citizens eligible to vote.}
    \item{votes:}{number of eligible citizens who did vote.}
    \item{invalid1, invalid2:}{number of invalid first and second votes (see
      details below).}
    \item{valid1, valid2:}{number of valid first and second votes.}
    \item{SPD1, SPD2:}{number of first and second votes for the Social Democrats.}
    \item{UNION1, UNION2:}{number of first and second votes for CDU/CSU,
      the conservative Christian Democrats.}
    \item{GRUENE1, GRUENE2:}{number of first and second votes for the
      Green Party.}
    \item{FDP1, FDP2:}{number of first and second votes for the Liberal Party.}
    \item{LINKE1, LINKE2:}{number of first and second votes for the Left
      Party (PDS in 2002).}
  }
  Missing values indicate that a party did not candidate in the
  corresponding constituency.
}
\arguments{
  \item{year}{numeric or character, year of the election.}
  \item{second}{logical, return second or first votes?}
  \item{percent}{logical, return percentages or absolute numbers?}
  \item{nazero}{logical, convert \code{NA}s to 0?}
  \item{state}{logical or character. If \code{TRUE} then only column \code{state}
    from the corresponding data frame is returned, and all other
    arguments are ignored. If character, then it is used as pattern to
    \code{\link{grep}} for the corresponding state(s), see examples.}
}
\details{
  \code{btw200x} are the original data sets.
  \code{bundestag()} is a helper function which extracts first
  or second votes, calculates percentages (number of votes for a party divided by
  number of valid votes), replaces missing values by zero, and converts
  the result from a data frame to a matrix.  By default
  it returns the percentage of second votes for each party, which
  determines the number of seats each party gets in parliament.
} 
\section{German Federal Elections}{ Half of the Members of the German
  Bundestag are elected directly from Germany's 299 constituencies, the
  other half on the parties' state lists.  Accordingly, each voter has
  two votes in the elections to the German Bundestag. The first vote,
  allowing voters to elect their local representatives to the Bundestag,
  decides which candidates are sent to Parliament from the
  constituencies.

  The second vote is cast for a party list. And it is this second vote
  that determines the relative strengths of the parties represented in
  the Bundestag.  At least 598 Members of the German Bundestag are
  elected in this way. In addition to this, there are certain
  circumstances in which some candidates win what
  are known as ``overhang mandates'' when the seats are being
  distributed.
}
\references{Homepage of the Bundestag:
  \url{http://www.bundestag.de}
}
\keyword{datasets}
\examples{
p02 <- bundestag(2002)
pairs(p02)
p05 <- bundestag(2005)
pairs(p05)
p09 <- bundestag(2009)
pairs(p09)

state <- bundestag(2002, state=TRUE)
table(state)

start.with.b <- bundestag(2002, state="^B")
table(start.with.b)

pairs(p09, col=2-(state=="Bayern"))
}
