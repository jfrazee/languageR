\name{nessdemog}
\alias{nessdemog}
\docType{data}
\title{Frequency spectrum for -ness in the context-governed BNC}
\description{
  Frequency (m) and frequency of frequency (Vm) for string types
  with the suffix \emph{-ness} in the demographic subcorpus of the 
  British National Corpus sampling spoken British English.
}
\usage{data(nessdemog)}
\format{
  A data frame with 15 observations on the following 2 variables.
  \describe{
    \item{\code{m}}{a numeric vector with word frequencies.}
    \item{\code{Vm}}{a numeric vector with the frequencies of word frequencies.}
  }
}
\source{
  The British National Corpus, see \url{http://www.natcorp.ox.ac.uk/}
}
\examples{
data(nessdemog)
library(zipfR)
nessdemog.spc = spc(m=nessdemog$m, Vm = nessdemog$Vm)
plot(nessdemog.spc)
}
\keyword{datasets}
