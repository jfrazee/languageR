\name{splitplot}
\alias{splitplot}
\docType{data}
\title{Simulated split-plot data set with subjects and items}
\description{
  Simulated lexical decision latencies with SOA as treatment, using
  a split-plot design with subjects and items, as available in
  Raaijmakers et al. (1999).
}
\usage{data(splitplot)}
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
data(splitplot)
items.split.plot.fnc(splitplot)
}}
\keyword{datasets}
