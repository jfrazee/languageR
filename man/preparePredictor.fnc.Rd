\name{preparePredictor.fnc}
\alias{preparePredictor.fnc}
\title{ determine X and Y values for a given (sub)plot}
\description{
  this function figures out the X and Y values for a given (sub)plot,
  including upper and lower 95\% HPD intervals
}
\usage{
preparePredictor.fnc(pred, model, m, ylabel, fun, val, xlabel, mcmc, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{ character string denoting predictor to be plotted on horizontal axis }
  \item{model}{ model fit by \code{lmer}}
  \item{m}{ matrix as produced by \code{makeDefaultMatrix.fnc} }
  \item{ylabel}{label for vertical axis (if other than name of dependent variable}
  \item{fun}{ character string denoting transformation function for dependent variable, currently only '"plogis"' or '"exp"' }
  \item{val}{ value of interacting variable }
  \item{xlabel}{ label for horizontal axis }
  \item{mcmc}{ MCMC matrix as produced by \code{mcmcsamp} or \code{pvals.fnc()$mcmc}}
  \item{\dots}{ further graphical parameters, currently not implemented }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  A data frame with values to be plotted, with columns
  \item{X}{values of predictor}
  \item{Y}{fitted values}
  \item{Type}{logical for whether predictor is factor}
  \item{Interaction}{logical for whether predictor is interacting predictor}
  \item{Levels}{for factors, the factor level names (only present for factors)}
}
\references{ ~put references to the literature/web site here ~ }
\author{R. H. Baayen}
\note{
 not intended for independent use
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
  \dontrun{not intended for independent use}
}
\keyword{ regression }
