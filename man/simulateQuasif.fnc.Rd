\name{simulateQuasif.fnc}
\alias{simulateQuasif.fnc}
\title{Simulate data for quasi-F analysis and compare models}
\description{
  This function creates a user-specified number of simulated  
  datasets, and compares mixed-effects models with quasi-F and 
	F1 and F2 analyses.
}
\usage{
simulateQuasif.fnc(dat, with = TRUE, nruns = 100, nsub = NA, nitem = NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame with a data set with as variables Subject, Item and
    SOA, as in the \code{quasif} data set.}
  \item{with}{Logical, if TRUE, an effect of SOA is built into the simulation.}
  \item{nruns}{Integer for the number of simulation runs.}
  \item{nsub}{Integer denoting the number of subjects.}
  \item{nitem}{Integer denoting the number of items.}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\details{
  Model parameters are estimated from the input data set.

  For each completed simulation run, a dot is added to the R console.
}
\value{
  A list with components
  \item{alpha05}{Description of 'comp1'}
  \item{alpha01}{proportion of runs in which predictors are significant at the
        05 significance level.}
  \item{res}{Data frame with simulation results.}
  \item{with}{Logical, TRUE if SOA effect is built into the simulations.}
}
\author{R. H. Baayen}
\seealso{See also \code{\link{subjects.quasif.fnc}}.}
\examples{\dontrun{
data(quasif)
library(lme4)

quasif.sim = simulateQuasif.fnc(quasif, nruns = 1000, with = TRUE) 
quasif.sim$alpha05
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
