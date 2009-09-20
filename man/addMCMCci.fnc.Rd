\name{addMCMCci.fnc}
\alias{addMCMCci.fnc}
\title{Add MCMC confidence intervals}
\description{
  Adds MCMC confidence intervals to the plot data frames for \code{plotLMER.fnc()}.
}
\usage{
addMCMCci.fnc(mcmcM, model, m, fun, pred, predname = NA, factor = FALSE)
}
\arguments{
  \item{mcmcM}{An MCMC matrix as obtained with \code{mcmcsamp()} or 
  \code{pvals.fnc()$mcmc}}
  \item{model}{A (generalized) linear mixed effects model obtained with
  \code{lmer()}}
  \item{m}{A design matrix constructed by \code{plotLMER.fnc()}}
  \item{fun}{A function for transforming the dependent variable}
  \item{pred}{The predictor to be plotted}
  \item{predname}{The name of the predictor}
  \item{factor}{Boolean for whether predictor is factor}
}
\details{
  Not intended for independent use.
}
\value{
  A data frame with components
  \item{x}{values for horizontal axis}
  \item{lower}{lower bound 95\% prediction interval}
  \item{upper}{upper bound of 95\% prediction interval}
}
\references{ ~put references to the literature/web site here ~ }
\author{R. H. Baayen}
\note{Under development}
\examples{
\dontrun{
  not intended for independent use.
}
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\keyword{ regression }

