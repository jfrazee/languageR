\name{spanishMeta}
\alias{spanishMeta}
\docType{data}
\title{Metadata for the spanish and spanishFunctionWords data sets}
\description{
  By-text metadata for the spanish and spanishFunctionWords data sets.
}
\usage{data(spanishMeta)}
\format{
  A data frame with 15 observations on the following 6 variables.
  \describe{
    \item{\code{Author}}{a factor with levels \code{C}, \code{M}, and \code{V}.}
    \item{\code{YearOfBirth}}{a numeric vector with year of birth of the author.}
    \item{\code{TextName}}{a factor with codes for the texts as levels (
    \code{X14458gll} ... \code{X14476gll}).}
    \item{\code{PubDate}}{a numeric vector with data of publication of the text.}
    \item{\code{Nwords}}{a numeric vector with text sizes in tokens.}
    \item{\code{FullName}}{a factor with author names: \code{Cela}, 
       \code{Mendoza} and \code{VargasLLosa}.}
  }
}
\references{
Spassova, M. S. (2006) \emph{Las marcas sintacticas de atribucion forense de
autoria de textos escritos en espanol}, Masters thesis, Institut Universitari
de Linguistica Aplicada, Universitat Pompeu Fabra, Barcelona. 
}
\examples{\dontrun{
data(spanishMeta)
}}
\keyword{datasets}
