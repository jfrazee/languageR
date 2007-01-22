\name{warlpiri}
\alias{warlpiri}
\docType{data}
\title{Ergative case marking in Warlpiri}
\description{
 This data set documents the use of ergative case marking in the narratives
 of native speakers of Lajamanu Warlpiri (8 children, 13 adults) 
 describing events in picture books.  
}
\usage{data(warlpiri)}
\format{
  A data frame with 347 observations on the following 9 variables.
  \describe{
    \item{\code{Speaker}}{a factor with speakers as levels.} 
    \item{\code{Sentence}}{a factor with sentence as levels.} 
    \item{\code{AgeGroup}}{a factor with levels \code{adult} and \code{child}.}
    \item{\code{CaseMarking}}{a factor with levels \code{ergative} and
       \code{other}.}
    \item{\code{WordOrder}}{a factor with levels \code{subInitial} (subject
       initial) and \code{subNotInitial} (subject not initial).}
    \item{\code{AnimacyOfSubject}}{a factor with levels \code{animate} 
       and \code{inanimate}.}
    \item{\code{OvertnessOfObject}}{a factor with levels \code{notOvert} 
       and \code{overt}.}
    \item{\code{AnimacyOfObject}}{a factor with levels \code{animate} 
       and \code{inanimate}.}
    \item{\code{Text}}{a factor with levels \code{texta}, 
       \code{textb} and \code{textc}.}
  }
}
\references{
O'Shannessy, C. (2006) \emph{Language contact and child bilingual acquisition:
Learning a mixed language and Warlpiri in northern Australia}, PhD Thesis,
University of Sydney, Australia.}
}
\examples{
data(warlpiri)
library(lme4, keep.source = FALSE)
warlpiri.lmer = lmer(CaseMarking ~ WordOrder * AgeGroup + 
  AnimacyOfSubject + (1|Text) + (1|Speaker), 
  family = "binomial", data = warlpiri, method = "Laplace")
warlpiri.lmer
}
\keyword{datasets}
