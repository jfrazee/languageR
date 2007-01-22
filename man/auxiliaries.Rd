\name{auxiliaries}
\alias{auxiliaries}
\docType{data}
\title{Auxiliaries for regular and irregular verbs in Dutch}
\description{
  For 285 regular and irregular Dutch verbs, the auxiliary for the present 
  and past perfect is listed together with the count
  of verbal synsets in WordNet.  Regular and irregular verbs are matched
  in the mean for lemma frequency.
}
\usage{data(auxiliaries)}
\format{
  A data frame with 285 observations on the following 4 variables.
  \describe{
    \item{\code{Verb}}{a factor with 285 monomorphemic Dutch verbs.}
    \item{\code{Aux}}{a factor with as levels the auxiliaries \code{hebben}, 
    \code{zijn} and \code{zijnheb} (for verbs allowing both auxiliaries).}
    \item{\code{VerbalSynsets}}{a numeric vector witth the number of 
          verbal synonym sets in WordNet in which the verb is listed.}
    \item{\code{Regularity}}{a factor with levels \code{irregular} and
          \code{regular}.}
  }
}
\references{
  Baayen, R. H. and Moscoso del Prado Martin, F. (2005)
  Semantic density and past-tense formation in three Germanic 
  languages, \emph{Language}, 81, 666-698.
}
\examples{
data(auxiliaries)
kruskal.test(auxiliaries$VerbalSynsets, auxiliaries$Aux)
}
\keyword{datasets}
