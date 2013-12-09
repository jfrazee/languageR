\name{sizeRatings}
\alias{sizeRatings}
\docType{data}
\title{Size ratings for 81 English concrete nouns}
\description{
  Subjective estimates of the size of the referents of 81 English 
  concrete nouns, collected from 38 subjects.
}
\usage{data(sizeRatings)}
\format{
  A data frame with 3078 observations on the following 7 variables.
  \describe{
    \item{\code{Rating}}{a numeric vector with subjective estimates of the
       size of the word's referent.}
    \item{\code{Subject}}{a factor with subjects as levels.}
    \item{\code{Word}}{a factor with words as levels.} 
    \item{\code{Class}}{a factor with levels \code{animal} and \code{plant}.}
    \item{\code{Naive}}{a factor with levels \code{naive} and \code{notNaive},
       coding whether the subject new about the purpose of the experiment.}
    \item{\code{Language}}{a factor with levels \code{English} and \code{notEnglish} coding whether the subject was a native speaker of English.}
    \item{\code{MeanFamiliarity}}{a numeric vector for the by-item mean
       familiarity ratings.}
  }
}
\details{
  Data collected with Jen Hay, University of Canterbury, Christchurch,
  New Zealand, 2004.
}
\examples{\dontrun{
data(sizeRatings)
require(lme4)
require(lmerTest)
require(optimx)
sizeRatings.lmer = lmer(Rating ~ Class * Naive + 
  MeanFamiliarity * Language + (1|Subject) + (1|Word), 
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")),
  data = sizeRatings)
summary(sizeRatings.lmer)
}}
\keyword{datasets}
