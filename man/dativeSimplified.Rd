\name{dativeSimplified}
\alias{dativeSimplified}
\docType{data}
\title{Dative Alternation - simplified data set}
\description{
Data describing the realization of the dative as NP or PP in the Switchboard
corpus and the Treebank Wall Street Journal collection. Simplified version
of the dative data set.
}
\usage{data(dativeSimplified)}
\format{
  A data frame with 903 observations on the following 5 variables.
  \describe{
    \item{\code{Verb}}{a factor with the verbs as levels.} 
    \item{\code{AnimacyOfRec}}{a factor with levels \code{animate} and
     \code{inanimate} for the animacy of the recipient.}
    \item{\code{LengthOfTheme}}{a numeric vector coding the number of words
      comprising the theme.}
    \item{\code{AnimacyOfTheme}}{a factor with levels \code{animate} and
      \code{inanimate} coding the animacy of the theme.}
    \item{\code{RealizationOfRec}}{a factor with levels \code{NP} and
      \code{PP} coding the realization of the dative.}
  }
}
\references{
Bresnan, J., Cueni, A., Nikitina, T. and Baayen, R. H. (2007) Predicting the
dative alternation, in Bouma, G. and Kraemer, I. and Zwarts, J.  (eds.),
\emph{Cognitive Foundations of Interpretation}, Royal Netherlands Academy
of Sciences, 33 pages, in press.
}
\examples{
data(dative)
}
\keyword{datasets}
