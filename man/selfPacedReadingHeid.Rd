\name{selfPacedReadingHeid}
\alias{selfPacedReadingHeid}
\docType{data}
\title{Self-paced reading latencies for Dutch neologisms}
\description{
  Self-paced reading latencies for Dutch neologisms ending in the
  suffix \emph{-heid}.
}
\usage{data(selfPacedReadingHeid)}
\format{
  A data frame with 1280 observations on the following 18 variables.
  \describe{
    \item{\code{Subject}}{a factor with subjects as levels.}
    \item{\code{Word}}{a factor with words as levels.}
    \item{\code{RT}}{a numeric vector with logarithmically transformed
      reading latencies.}
    \item{\code{RootFrequency}}{a numeric vector for the logarithmically
      transformed frequency of the lowest-level base of the neologism
      (e.g., \emph{lob} in \emph{[[[lob]+ig]+heid]}.}
    \item{\code{Condition}}{a factor for the priming conditions 
      with levels \code{baseheid} (neologism is preceded 40 trials back
      by its base word) and \code{heidheid} (the neologism is preceded
      40 trials back by itself).}
    \item{\code{Rating}}{a numeric vector for the word's subjective frequency
      estimate.}
    \item{\code{Frequency}}{a numeric vector for the neologism's frequency
      (all zero).}
    \item{\code{BaseFrequency}}{a numeric vector for the base adjective
      underlying the neologism (e.g., \emph{lobbig} in 
      \emph{[[[lob]+ig]+heid]}).}
    \item{\code{LengthInLetters}}{a numeric vector coding word length in
      letters.}
    \item{\code{FamilySize}}{a numeric vector for the logaritmically
      transformed count of a word's morphological family members.}
    \item{\code{NumberOfSynsets}}{a numeric vector for the count of synonym
      sets in WordNet in which the word is listed.}
    \item{\code{RT4WordsBack}}{a numeric vector for the log-transformed 
      reading latencies four trials back.}
    \item{\code{RT3WordsBack}}{a numeric vector for the log-transformed 
      reading latencies three trials back.}
    \item{\code{RT2WordsBack}}{a numeric vector for the log-transformed 
      reading latencies two trials back.}
    \item{\code{RT1WordBack}}{a numeric vector for the log-transformed 
      reading latencies one trial back.}
    \item{\code{RT1WordLater}}{a numeric vector for the log-transformed 
      reading latencies one trial later.}
    \item{\code{RT2WordsLater}}{a numeric vector for the log-transformed 
      reading latencies two trials later.}
    \item{\code{RTtoPrime}}{a numeric vector for the log-transformed reading
     latency for the prime.}
  }
}
\references{
De Vaan, L., Schreuder, R. and Baayen, R. H. (2007) Regular morphologically
complex neologisms leave detectable traces in the mental lexicon, \emph{The
Mental Lexicon}, 2, in press.
}
\examples{\dontrun{
data(selfPacedReadingHeid)

# data validation
plot(sort(selfPacedReadingHeid$RT))   
selfPacedReadingHeid = selfPacedReadingHeid[selfPacedReadingHeid$RT > 5 & 
  selfPacedReadingHeid$RT < 7.2,]

# fitting a mixed-effects model

library(lme4, keep.source = FALSE)
x = selfPacedReadingHeid[,12:15]
x.pr = prcomp(x, center = TRUE, scale = TRUE)
selfPacedReadingHeid$PC1 = x.pr$x[,1]

selfPacedReadingHeid.lmer = lmer(RT ~ RTtoPrime + LengthInLetters + 
  PC1 * Condition + (1|Subject) + (1|Word), 
  data = selfPacedReadingHeid)  
pvals.fnc(selfPacedReadingHeid.lmer)$summary

# model criticism

selfPacedReadingHeid.lmerA = lmer(RT ~ RTtoPrime + LengthInLetters + 
PC1 * Condition + (1|Subject) + (1|Word), data = 
selfPacedReadingHeid[abs(scale(resid(selfPacedReadingHeid.lmer))) < 2.5, ])

qqnorm(resid(selfPacedReadingHeid.lmerA))

pvals.fnc(selfPacedReadingHeid.lmerA)$summary
}}
\keyword{datasets}
