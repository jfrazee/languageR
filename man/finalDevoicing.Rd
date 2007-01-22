\name{finalDevoicing}
\alias{finalDevoicing}
\docType{data}
\title{Final Devoicing in Dutch}
\description{
Phonological specifications for onset, nucleus and offset for 1697 Dutch
monomorphemic words with a final obstruent.  These final obstruents may exhibit
a voicing alternation that is traditionally described as syllable-final
devoicing: underlying /d/ in /hond/ becomes a /t/ when syllable-final ([hOnt])
and remains a /d/ otherwise ([hOn-den]). 
}
\usage{data(finalDevoicing)}
\format{
  A data frame with 1697 observations on the following 9 variables.
  \describe{
    \item{\code{Word}}{a factor with the words as levels.} 
    \item{\code{Onset1Type}}{a factor for the first consonant in the onset, with levels \code{None}, 
          \code{Obstruent} and \code{Sonorant}.}
    \item{\code{Onset2Type}}{a factor for the second consonant in the onset, with levels \code{None}, 
          \code{Obstruent} and \code{Sonorant}.}
    \item{\code{VowelType}}{a factor describing the vowel with levels \code{iuy}, \code{long} and 
       \code{short}.}
    \item{\code{ConsonantType}}{a factor for the first consonant in the offset, with levels \code{None},
          \code{Obstruent} and \code{Sonorant}.}
    \item{\code{Obstruent}}{a factor describing place and manner of articulation of the final obstruent, 
          with levels \code{F} (/f,v/), \code{P} (/p,b/), \code{S} (/s,z/), \code{T} (/t,d/) and
          \code{X} (/x,g/).}
    \item{\code{Nsyll}}{a numeric vector for the number of syllables in the word.}
    \item{\code{Stress}}{a factor with levels \code{A} (antepenult), \code{F} (final) and
          \code{P} (penult).}
    \item{\code{Voice}}{a factor with levels \code{voiced} and \code{voiceless}.}
  }
}
\references{
Ernestus, M. and Baayen, R. H. (2003) Predicting the unpredictable:
Interpreting neutralized segments in Dutch, \emph{Language}, 79, 5-38.
}
\examples{
data(finalDevoicing)
library(rpart)

# ---- CART tree 

finalDevoicing.rp = rpart(Voice ~ ., data = finalDevoicing[ , -1])
plotcp(finalDevoicing.rp)
finalDevoicing.pruned = prune(finalDevoicing.rp, cp = 0.021)
plot(finalDevoicing.pruned, margin = 0.1, compress = TRUE)
text(finalDevoicing.pruned, use.n = TRUE, pretty = 0, cex=0.8)

# ---- logistic regression 

library(Design)

finalDevoicing.dd = datadist(finalDevoicing)
options(datadist='finalDevoicing.dd')

finalDevoicing.lrm = lrm(Voice ~ VowelType + ConsonantType + Obstruent + 
Nsyll + Stress + Onset1Type + Onset2Type, data = finalDevoicing)
anova(finalDevoicing.lrm)

# ---- model simplification

fastbw(finalDevoicing.lrm)

finalDevoicing.lrm = lrm(Voice ~ VowelType + ConsonantType + 
Obstruent + Nsyll, data = finalDevoicing, x = TRUE, y = TRUE)

par(mfrow=c(2,2))
plot(finalDevoicing.lrm, fun = plogis, ylim = c(0, 1), 
ylab = "p(voiceless)")
par(mfrow=c(1,1))

# ---- model validation

validate(finalDevoicing.lrm, B = 200)
}
\keyword{datasets}
