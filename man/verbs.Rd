\name{verbs}
\alias{verbs}
\docType{data}
\title{Dative Alternation - simplified data set}
\description{
A simplified version of the \code{dative} data set, used
for expository purposes only.
}
\usage{data(verbs)}
\format{
  A data frame with 903 observations on the following 5 variables.
  \describe{
    \item{\code{RealizationOfRec}}{a factor with levels \code{NP} 
          and \code{PP}.}
    \item{\code{Verb}}{a factor with the verbs as levels.} 
    \item{\code{AnimacyOfRec}}{a factor with levels \code{animate} and
          \code{inanimate}.}
    \item{\code{AnimacyOfTheme}}{a factor with levels \code{animate} and
          \code{inanimate}.}
    \item{\code{LengthOfTheme}}{a numeric vector coding the length in words
          of the theme.}
  }
}
\references{
Bresnan, J., Cueni, A., Nikitina, T. and Baayen, R. H. (2007) Predicting the
dative alternation, in Bouma, G. and Kraemer, I. and Zwarts, J.  (eds.),
\emph{Cognitive Foundations of Interpretation}, Royal Netherlands Academy
of Sciences, 33 pages, in press.
}
\examples{
data(verbs)
head(verbs)
xtabs( ~ RealizationOfRec + AnimacyOfRec, data = verbs)
barplot(xtabs( ~ RealizationOfRec + AnimacyOfRec, data = verbs),beside=TRUE)
}
\keyword{datasets}
