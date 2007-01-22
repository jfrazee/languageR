\name{durationsOnt}
\alias{durationsOnt}
\docType{data}
\title{Durational measurements on the Dutch prefix ont-}
\description{
Durational measurements on the Dutch prefix \emph{ont-} 
in the Spoken Dutch Corpus. 
}
\usage{data(durationsOnt)}
\format{
  A data frame with 102 observations on the following 11 variables.
  \describe{
    \item{\code{Word}}{a factor with the words as levels.}
    \item{\code{Frequency}}{a numeric vector with the word's logarithmically
       transformed frequency in the Spoken Dutch Corpus.}
    \item{\code{Speaker}}{a factor with speakers as levels.}
    \item{\code{Sex}}{a factor with levels \code{female} and \code{male}.}
    \item{\code{YearOfBirth}}{a numeric vector coding year of birth of the
       speaker - 1900.}
    \item{\code{DurationOfPrefix}}{a numeric vector for the duration of ont-
       in seconds}
    \item{\code{DurationPrefixVowel}}{a numeric vector for the duration of the
       vowel in the prefix in seconds.}
    \item{\code{DurationPrefixNasal}}{a numeric vector for the duration of the
       nasal in the prefix in seconds.}
    \item{\code{DurationPrefixPlosive}}{a numeric vector for the duration of the   plosive in the prefix in seconds.}
    \item{\code{NumberOfSegmentsOnset}}{a numeric vector for the number of segments in the onset of the stem.}
    \item{\code{PlosivePresent}}{a factor with levels \code{no} and \code{yes} 
       for whether the plosive is realized in the signal.}
    \item{\code{SpeechRate}}{a numeric vector coding speech rate in number
       of syllables per second.}
  }
}
\references{
Pluymaekers, M., Ernestus, M. and Baayen, R. H. (2005) Frequency and acoustic
length: the case of derivational affixes in Dutch, \emph{Journal of the
Acoustical Society of America}, 118, 2561-2569.
}
\examples{
data(durationsOnt)

###### modeling the duration of the prefix

prefix.lm = lm(DurationOfPrefix ~ (YearOfBirth + SpeechRate) * Frequency, 
   data = durationsOnt)
summary(prefix.lm)

# ---- model criticism

plot(prefix.lm)
outliers = c(36, 35, 17, 72)
prefix.lm = lm(DurationOfPrefix ~ (YearOfBirth + SpeechRate) * Frequency, 
   data = durationsOnt[-outliers,])
summary(prefix.lm)

###### modeling the presence of the /t/ 

library(Design)
durationsOnt.dd = datadist(durationsOnt)
options(datadist = 'durationsOnt.dd')

plosive.lrm = lrm(PlosivePresent ~ SpeechRate + YearOfBirth, 
   data = durationsOnt, x = TRUE, y = TRUE)
plosive.lrm
validate(plosive.lrm, bw = TRUE, B = 200)

###### modeling the duration of the /n/

nasal.lm = lm(DurationPrefixNasal ~ PlosivePresent + Frequency + 
   YearOfBirth, data = durationsOnt)
summary(nasal.lm)

# ---- model criticism

plot(nasal.lm)
outliers = c(71, 28, 62, 33)
nasal.lm = lm(DurationPrefixNasal ~ PlosivePresent + Frequency + 
   YearOfBirth, data = durationsOnt[-outliers,])
summary(nasal.lm)

}
\keyword{datasets}
