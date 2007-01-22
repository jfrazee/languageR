\name{latinsquare}
\alias{latinsquare}
\docType{data}
\title{Simulated Latin Square data set with subjects and items}
\description{
  Simulated lexical decision latencies with SOA as treatment, using
  a Latin Square design with subjects and items, as available in
  Raaijmakers et al. (1999).
}
\usage{data(latinsquare)}
\format{
  A data frame with 144 observations on the following 6 variables.
  \describe{
    \item{\code{Group}}{a factor with levels \code{G1}, \code{G2} and 
        \code{G3}, for groups of subjects}
    \item{\code{Subject}}{a factor with subjects labelled
    \code{S1}, ... \code{S12}.}
    \item{\code{Word}}{a factor with words labelled \code{W1} ... \code{W12}.}
    \item{\code{RT}}{a numeric vector for reaction times.}
    \item{\code{SOA}}{a factor with levels \code{long}, \code{medium}, 
        and \code{short}.}
    \item{\code{List}}{a factor with levels \code{L1}, \code{L2}, and \code{L3}
        for lists of words.}
  }
}
\source{
Raaijmakers, J.G.W., Schrijnemakers, J.M.C. & Gremmen, F. (1999)
How to deal with "The language as fixed effect fallacy": 
common misconceptions and alternative solutions, 
\emph{Journal of Memory and Language}, 41, 416-426.
}
\examples{
data(latinsquare)
\dontrun{
library(lme4)
latinsquare.with = 
   simulate.latinsquare.fnc(latinsquare, nruns = 1000, with = TRUE) 
latinsquare.without = 
   simulate.latinsquare.fnc(latinsquare, nruns = 1000, with = FALSE)
latinsquare.with$alpha05
latinsquare.without$alpha05
}
}
\keyword{datasets}
