\name{writtenVariationLijk}
\alias{writtenVariationLijk}
\docType{data}
\title{Variation in written Dutch in the use of the suffix -lijk}
\description{
  This dataset documents variation in the use of the 80 most frequent
  words ending in the suffix \emph{-lijk} in written Dutch.
}
\usage{data(writtenVariationLijk)}
\format{
  A data frame with 560 observations on the following 5 variables.
  \describe{
    \item{\code{Corpus}}{a factor with as levels the sampled newspapers: 
       \code{belang} (Het Belang van Limburg), 
       \code{gazet} (De Gazet van Antwerpen),
       \code{laatnieu} (Het Laatste Nieuws),
       \code{limburg} (De Limburger),
       \code{nrc} (NRC Handelsblad),
       \code{stand} (De Standaard), and
       \code{tele} (De Telegraaf).}
    \item{\code{Word}}{a factor with the 80 most frequent words ending in
          \emph{-lijk}.}
    \item{\code{Count}}{a numeric vector with token counts in the CONDIV
        corpus.}
    \item{\code{Country}}{a factor with levels \code{Flanders} and 
       \code{Netherlands}.}
    \item{\code{Register}}{a factor with levels \code{National}, 
       \code{Quality} and \code{Regional} coding the type of newspaper.}
  }
}
\references{
Keune, K., Ernestus, M., Van Hout, R. and Baayen, R.H. (2005) Social,
geographical, and register variation in Dutch: From written 'mogelijk' to
spoken 'mok', \emph{Corpus Linguistics and Linguistic Theory}, 1, 183-223.
}
\examples{\dontrun{
data(writtenVariationLijk)

require(lme4)
require(lmerTest)
require(lme4)

writtenVariationLijk.lmer = glmer(Count ~ Country * Register + (1|Word), 
  control=glmerControl(optimizer="optimx",optCtrl=list(method="nlminb")),
  data = writtenVariationLijk, family = "poisson")

writtenVariationLijk.lmerA = glmer(Count ~ Country * Register + (Country|Word), 
  control=glmerControl(optimizer="optimx",optCtrl=list(method="nlminb")),
  data = writtenVariationLijk, family = "poisson")

anova(writtenVariationLijk.lmer, writtenVariationLijk.lmerA)

summary(writtenVariationLijk.lmerA)
}}
\keyword{datasets}
