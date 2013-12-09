\name{primingHeid}
\alias{primingHeid}
\docType{data}
\title{Primed lexical decision latencies for neologisms ending in -heid}
\description{
  Primed lexical decision latencies for Dutch neologisms ending in the
  suffix \emph{-heid}.
}
\usage{data(primingHeid)}
\format{
  A data frame with 832 observations on the following 13 variables.
  \describe{
    \item{\code{Subject}}{a factor with subjects as levels.}
    \item{\code{Word}}{a factor with words as levels.}
    \item{\code{Trial}}{a numeric vector for the rank of the trial in its
       experimental list.}
    \item{\code{RT}}{a numeric vector with log-transformed lexical decision
       latencies.}
    \item{\code{Condition}}{a factor coding the priming treatmen,
       with levels \code{baseheid} (prime is the base word) and
       \code{heid} (the prime is the neologism)}
    \item{\code{Rating}}{a numeric vector for subjective frequency estimates.}
    \item{\code{Frequency}}{a numeric vector for 
       log-transformed frequencies of the whole word.}
    \item{\code{BaseFrequency}}{a numeric vector for the log-transformed 
       frequencies of the base word.}
    \item{\code{LengthInLetters}}{a numeric vector coding orthographic length
       in letters.}
    \item{\code{FamilySize}}{a numeric vector for the log-transformed
       count of the word's morphological family.}
    \item{\code{NumberOfSynsets}}{a numeric vector for the number of synonym
       sets in WordNet in which the base is listed.}
    \item{\code{ResponseToPrime}}{a factor with levels \code{correct} and
       \code{incorrect} for the response to the prime.}
    \item{\code{RTtoPrime}}{a numeric vector for the log-transformed 
       reaction time to the prime.}
  }
}
\references{
De Vaan, L., Schreuder, R. and Baayen, R. H. (2007) Regular morphologically
complex neologisms leave detectable traces in the mental lexicon, \emph{The
Mental Lexicon}, 2, in press.
}
\examples{\dontrun{
data(primingHeid)

require(lme4)
require(lmerTest)
require(optimx)

primingHeid.lmer = lmer(RT ~ RTtoPrime * ResponseToPrime + Condition +
  (1|Subject) + (1|Word), 
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")),
  data = primingHeid)
summary(primingHeid.lmer)

# model criticism

primingHeid.lmer = lmer(RT ~ RTtoPrime * ResponseToPrime + Condition +
  (1|Subject) + (1|Word), 
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")),
  data = primingHeid[abs(scale(resid(primingHeid.lmer)))<2.5,])
summary(primingHeid.lmer)
}}
\keyword{datasets}
