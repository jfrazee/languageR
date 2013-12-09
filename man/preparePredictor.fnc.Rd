\name{preparePredictor.fnc}
\alias{preparePredictor.fnc}
\title{ determine X and Y values for a given (sub)plot}
\description{
  this function figures out the X and Y values for a given (sub)plot,
  including upper and lower 95\% HPD intervals
}
\usage{
preparePredictor.fnc(pred, model, m, ylabel, fun, val, xlabel, ranefs, ...)
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
  \item{ranefs}{ a three-element list \code{Group, Level, Predictor}, specifying a random-effect Group (e.g. \code{Subject}), a level (e.g., \code{S10}) and a value (e.g., \code{LogFrequency}) for which partial effects have to be calibrated; implemented only for \code{mcmcMat=NA}.}
  \item{\dots}{ further graphical parameters, currently not implemented }
}
%\details{
%  Not intended for independent use.
%}
\value{
  A data frame with values to be plotted, with columns
  \item{X}{values of predictor}
  \item{Y}{fitted values}
  \item{Type}{logical for whether predictor is factor}
  \item{Interaction}{logical for whether predictor is interacting predictor}
  \item{Levels}{for factors, the factor level names (only present for factors)}
}
%\references{None}
\author{R. H. Baayen}
\note{
 not intended for independent use
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
  \dontrun{Not intended for independent use.}
}
\keyword{ regression }
