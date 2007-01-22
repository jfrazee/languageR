\name{durationsGe}
\alias{durationsGe}
\docType{data}
\title{Durational measurements on the Dutch prefix ge-}
\description{
Durational measurements on the Dutch prefix \emph{ge-} in the 
Spoken Dutch Corpus.
}
\usage{data(durationsGe)}
\format{
  A data frame with 428 observations on the following 8 variables.
  \describe{
    \item{\code{Word}}{a factor with the words as levels.} 
    \item{\code{Frequency}}{a numeric vector with the word's 
       absolute frequency in the Spoken Dutch Corpus.}
    \item{\code{Speaker}}{a factor with the speakers as levels.}
    \item{\code{Sex}}{a factor with levels \code{female} and \code{male}, 
       this information is missing for one speaker.}
    \item{\code{YearOfBirth}}{a numeric vector with years of birth.}
    \item{\code{DurationOfPrefix}}{a numeric vector with the duration of the
         prefix -ont in seconds.}
    \item{\code{SpeechRate}}{a numeric vector coding speech rate in number
       of syllables per second.}
    \item{\code{NumberSegmentsOnset}}{a numeric vector for the number of 
          segments in the onset of the stem.}
  }
}
\references{
Pluymaekers, M., Ernestus, M. and Baayen, R. H. (2005) Frequency and acoustic
length: the case of derivational affixes in Dutch, \emph{Journal of the
Acoustical Society of America}, 118, 2561-2569.
}
\examples{
data(durationsGe)
durationsGe$Frequency = log(durationsGe$Frequency + 1)
durationsGe$YearOfBirth = durationsGe$YearOfBirth - 1900

durationsGe.lm = lm(DurationOfPrefix ~ Frequency+SpeechRate, data = durationsGe)
summary(durationsGe.lm)

# ---- model criticism

plot(durationsGe.lm)
outliers = c(271, 392, 256, 413, 118, 256)
durationsGe.lm = lm(DurationOfPrefix ~ Frequency + SpeechRate, 
  data = durationsGe[-outliers, ])
summary(durationsGe.lm)

}
\keyword{datasets}
