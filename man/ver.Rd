\name{ver}
\alias{ver}
\docType{data}
\title{The Dutch prefix ver-: semantic transparency and frequency}
\description{
  Semantic transparency (dichotomous) and frequency for 985 words
  with the Dutch prefix \emph{ver-}.
}
\usage{data(ver)}
\format{
  A data frame with 985 observations on the following 2 variables.
  \describe{
    \item{\code{Frequency}}{a numeric vector for the words' frequency.}
    \item{\code{SemanticClass}}{a factor with levels \code{opaque} and
      \code{transparent} coding semantic transparency.}
  }
}
\references{
Baayen, R. H. and Lieber, R. (1997) Word Frequency Distributions and Lexical
Semantics, \emph{Computers and the {H}umanities}, 30, 281-291.
}
\examples{\dontrun{
data(ver)
ver$Frequency = log(ver$Frequency)

plot(density(ver$Frequency))

# plot separate densities for opaque and transparent words

ver.transp = ver[ver$SemanticClass == "transparent",]$Frequency
ver.opaque = ver[ver$SemanticClass == "opaque", ]$Frequency

ver.transp.d = density(ver.transp)
ver.opaque.d = density(ver.opaque)
xlimit = range(ver.transp.d$x, ver.opaque.d$x)
ylimit = range(ver.transp.d$y, ver.opaque.d$y)
plot(ver.transp.d, lty = 1, col = "black", 
  xlab = "frequency", ylab = "density", 
  xlim = xlimit, ylim = ylimit, main = "")
lines(ver.opaque.d, col = "darkgrey")
legend(6,0.25, lty=rep(1,2), col=c("black", "darkgrey"), 
legend=c("transparent", "opaque"))

# test whether the difference is significant

ks.test(jitter(ver.transp), jitter(ver.opaque))
}}
\keyword{datasets}
