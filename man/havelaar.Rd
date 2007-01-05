\name{havelaar}
\alias{havelaar}
\docType{data}
\title{The determiner 'het' in the Dutch novel Max Havelaar}
\description{
  The frequency of the determiner 'het' in the Dutch novel
  'Max Havelaar' by Multatuli (Eduard Douwes Dekker), in 99
  consecutive text fragments of 1000 tokens each.
}
\usage{data(havelaar)}
\format{
  A data frame with 99 observations on the following 2 variables.
  \describe{
    \item{\code{Chunk}}{a numeric vector with the indices of the text 
    fragments.}
    \item{\code{Frequency}}{a numeric vector with the frequencies of the 
    determiner 'het' in the text fragments.}
  }
}
\source{
The text of Max Havelaar was obtained from the Project Gutenberg at
at \url{http://www.gutenberg.org/wiki/Main_Page}
}
\examples{\dontrun{
data(havelaar)

n = 1000                          # token size of text fragments
p = mean(havelaar$Frequency / n)  # relative frequencies

plot(qbinom(ppoints(99), n, p), sort(havelaar$Frequency),
   xlab = paste("quantiles of (", n, ",", round(p, 4), 
   ")-binomial", sep=""), ylab = "frequencies")


lambda = mean(havelaar$Frequency)
ks.test(havelaar$Frequency, "ppois", lambda)
ks.test(jitter(havelaar$Frequency), "ppois", lambda)


}}
\keyword{datasets}
