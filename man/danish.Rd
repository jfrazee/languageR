\name{danish}
\alias{danish}
\docType{data}
\title{Danish auditory lexical decision}
\description{
Auditory lexical decision latencies for Danish complex words.
}
\usage{data(danish)}
\format{
  A data frame with 3326 observations on the following 16 variables.
  \describe{
    \item{\code{Subject}}{a random-effect factor coding participants in the
       experiment.}
    \item{\code{Word}}{a random-effect factor coding the words for which 
       auditory lexical decision responses were elicited.}
    \item{\code{Affix}}{a random-effect factor coding the affixes 
       in the words.}
    \item{\code{LogRT}}{the dependent variable, log response latency.} 
    \item{\code{PC1}}{first principal component orthogonalizing the four 
      response latencies preceding the current trial in the experiment.}
    \item{\code{PC2}}{second principal component orthogonalizing the four 
      response latencies preceding the current trial in the experiment.}
    \item{\code{PrevError}}{factor with levels \code{CORRECT} and \code{ERROR} 
      coding whether the preceding trial elicited a correct lexical decision.} 
    \item{\code{Rank}}{the trial number in the experiment.} 
    \item{\code{Sex}}{factor coding the sex of the participant, 
      with levels \code{F} (female) and \code{M} (male).} 
    \item{\code{LogWordFreq}}{log-transformed word frequency.} 
    \item{\code{LogAffixFreq}}{log-transformed affix frequency.} 
    \item{\code{ResidFamSize}}{residualized morphological family size (taking out \code{LogWordFreq} and \code{LogAffixFreq}).} 
    \item{\code{ResidSemRating}}{residualized semantic rating (taking out
      morphological family size).} 
    \item{\code{LogCUP}}{log-transformed complex uniqueness point (CUP).} 
    \item{\code{LogUP}}{log-transformed uniqueness point (UP).}
    \item{\code{LogCUPtoEnd}}{log of the distance (in msec) from the CUP to
      the end of the word.} 
  }
}
\references{
L. Balling and R. H. Baayen (2008). Morphological effects in auditory word 
recognition: Evidence from Danish. Submitted to Language and Cognitive Processes.
}
\examples{
\dontrun{
data(danish)
require(lme4, quietly=TRUE)

# ---- mixed-effects regression with three random intercepts

danish.lmer = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
(1|Subject) + (1|Word) + (1|Affix), data = danish)
 
danish.lmerA = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
(1|Subject) + (1|Word) + (1|Affix), data = danish,
subset=abs(scale(resid(danish.lmer)))<2.5)

mcmc = pvals.fnc(danish.lmerA, nsim=10000, withMCMC=TRUE)
mcmc$fixed[,1:5]
}
}
\keyword{datasets}
