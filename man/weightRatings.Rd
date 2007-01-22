\name{weightRatings}
\alias{weightRatings}
\docType{data}
\title{Subjective estimates of the weight of the referents of 81 English nouns}
\description{
  Subjective estimates on a seven-point scale of the weight of the
  referents of 81 English nouns.
}
\usage{data(weightRatings)}
\format{
  A data frame with 1620 observations on the following 7 variables.
  \describe{
    \item{\code{Subject}}{a factor with subjects as levels.} 
    \item{\code{Rating}}{a numeric vector.}
    \item{\code{Trial}}{a numeric vector with the weight ratings.}
    \item{\code{Sex}}{a factor with levels \code{F} and \code{M}.}
    \item{\code{Word}}{a factor with words as levels.} 
    \item{\code{Frequency}}{a numeric vector with log-transformed
      lemma frequencies as available in the CELEX lexical database.}
    \item{\code{Class}}{a factor with levels \code{animal} and \code{plant}.}
  }
}
\references{
Data collected with Jen Hay, University of Canterbury, Christchurch,
New Zealand, 2004.
}
\examples{
data(weightRatings)
xylowess.fnc(Rating ~ Frequency | Subject, data = weightRatings, 
  xlab = "log Frequency", ylab = "Weight Rating")
}
\keyword{datasets}
