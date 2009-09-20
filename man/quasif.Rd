\name{quasif}
\alias{quasif}
\docType{data}
\title{Simulated data set with subjects and items requiring quasi-F ratios}
\description{
  Simulated lexical decision latencies with SOA as treatment, traditionally
	requiring an analysis using quasi-F ratios, as available in
  Raaijmakers et al. (1999).
}
\usage{data(quasif)}
\format{
  A data frame with 64 observations on the following 4 variables.
  \describe{
    \item{\code{Subject}}{a factor coding subjects.}
    \item{\code{RT}}{a numeric vector for simulated reaction times 
            in lexical decision.}
    \item{\code{Item}}{a factor coding items.}
    \item{\code{SOA}}{a factor coding SOA treatment with levels \code{long}
            and \code{short}.} 
  }
}
\source{
Raaijmakers, J.G.W., Schrijnemakers, J.M.C. & Gremmen, F. (1999)
How to deal with "The language as fixed effect fallacy": 
common misconceptions and alternative solutions, 
\emph{Journal of Memory and Language}, 41, 416-426.
}
\examples{\dontrun{
data(quasif)
items.quasif.fnc(quasif)
}}
\keyword{datasets}
