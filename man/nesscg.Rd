\name{nesscg}
\alias{nesscg}
\docType{data}
\title{Frequency spectrum for -ness in the demographic BNC}
\description{
  Frequency (m) and frequency of frequency (Vm) for string types
  with the suffix \emph{-ness} in the context-governed subcorpus of the 
  British National Corpus sampling spoken British English.
}
\usage{data(nesscg)}
\format{
  A data frame with 37 observations on the following 2 variables.
  \describe{
    \item{\code{m}}{a numeric vector with word frequencies.}
    \item{\code{Vm}}{a numeric vector with the frequencies of word frequencies.}
  }
}
\source{
  The British National Corpus, see \url{http://www.natcorp.ox.ac.uk/}
}
\examples{\dontrun{
data(nesscg)
library(zipfR)
nesscg.spc = spc(m=nesscg$m, Vm = nesscg$Vm)
plot(nesscg.spc)
}}
\keyword{datasets}
