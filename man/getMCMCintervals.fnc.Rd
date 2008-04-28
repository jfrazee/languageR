\name{getMCMCintervals.fnc}
\alias{getMCMCintervals.fnc}
\title{calculate HPD prediction intervals}
\description{
  calculate HPD 95\% prediction intervals
}
\usage{
getMCMCintervals.fnc(fixf, mcmcMatrix, m)
}
\arguments{
  \item{fixf}{ vector of fixed effects coefficients (\code{fixef(model.lmer)}) }
  \item{mcmcMatrix}{ MCMC matrix obtained with \code{mcmcsamp} or 
    \code{pvals.fnc}}
  \item{m}{ model matrix }
}
\details{
  not intended for independent use
}
\value{
  A matrix with columns '"lower"' and
  '"upper"' and rows corresponding to the values of the predictor to
  be plotted on the X-axis.
}
\references{languageR}
\author{R. H. Baayen}
\note{
not intended for independent use
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
  \dontrun{not intended for independent use}
}
\keyword{ regression }
