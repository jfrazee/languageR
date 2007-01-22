\name{heid}
\alias{heid}
\docType{data}
\title{Lexical decision latencies for words ending in -heid}
\description{
  A simplified version of the \code{primingHeid} dataset.
}
\usage{data(heid)}
\format{
  A data frame with 832 observations on the following 4 variables.
  \describe{
    \item{\code{Subject}}{a factor with subjects as levels.} 
    \item{\code{Word}}{a factor with words as levels.} 
    \item{\code{RT}}{a numeric vector with logarithmically transformed
       reaction times in visual lexical decision.}
    \item{\code{BaseFrequency}}{a numeric vector with the logarithmically
       transformed frequency of the base adjective of the word with the
       suffix \emph{-heid}.}
  }
}
\references{
De Vaan, L., Schreuder, R. and Baayen, R. H. (2007) Regular morphologically
complex neologisms leave detectable traces in the mental lexicon, \emph{The
Mental Lexicon}, 2, in press.
}
\examples{
data(heid)
heid = aggregate(heid$RT, list(heid$Word, heid$BaseFrequency), mean)
colnames(heid) = c("Word", "BaseFrequency", "MeanRT")
}
\keyword{datasets}
