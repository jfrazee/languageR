\name{twente}
\alias{twente}
\docType{data}
\title{Frequency spectrum for the Twente News Corpus}
\description{
  Frequency (m) and frequency of frequency (Vm) for string types
  in the Twente News Corpus.
}
\usage{data(twente)}
\format{
  A data frame with 4639 observations on the following 2 variables.
  \describe{
    \item{\code{m}}{a numeric vector with word frequencies.}
    \item{\code{Vm}}{a numeric vector with the frequencies of word frequencies.}
  }
}
\source{
  Twente News Corpus.
}
\examples{\dontrun{
data(twente)
library(zipfR)
twente.spc = spc(m=twente$m, Vm = twente$Vm)
plot(twente.spc)
}}
\keyword{datasets}
