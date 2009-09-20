\name{nessw}
\alias{nessw}
\docType{data}
\title{Frequency spectrum for -ness in the written BNC}
\description{
  Frequency (m) and frequency of frequency (Vm) for string types
  with the suffix \emph{-ness} in the subcorpus of the British
  National Corpus sampling written British English.
}
\usage{data(nessw)}
\format{
  A data frame with 189 observations on the following 2 variables.
  \describe{
    \item{\code{m}}{a numeric vector with word frequencies.}
    \item{\code{Vm}}{a numeric vector with the frequencies of word frequencies.}
  }
}
\source{
  The British National Corpus, see \url{http://www.natcorp.ox.ac.uk/}
}
\examples{\dontrun{
data(nessw)
library(zipfR)
nessw.spc = spc(m=nessw$m, Vm = nessw$Vm)
plot(nessw.spc)
}}
\keyword{datasets}
