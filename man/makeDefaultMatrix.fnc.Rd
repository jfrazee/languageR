\name{makeDefaultMatrix.fnc}
\alias{makeDefaultMatrix.fnc}
\title{ Create model matrix with main effects only }
\description{
  Creates a model matrix with main effects only
}
\usage{
makeDefaultMatrix.fnc(model, n = 100, conditioningPred = "", conditioningValue = NULL, control = NA)
}
\arguments{
  \item{model}{A model fit with \code{lmer}}
  \item{n}{integer specifying number of points to be plotted on X-axis}
  \item{conditioningPred}{name of predictor entering into interaction}
  \item{conditioningValue}{vector of values (numeric or factor level names) to be shown for interaction }
  \item{control}{a two-element list (predictor, value) specifying an additional predictor to be fixed to the given value in a partial effect plot. May be useful for hand-made plots for three-way interactions.}
}
\details{
  not intended for independent use
}
\value{
  a (model) matrix
}
%\references{ ~put references to the literature/web site here ~ }
\author{R. H. Baayen}
\note{ 
  not intended for independent use
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
  \dontrun{not intended for independent use}
}
\keyword{ regression }
