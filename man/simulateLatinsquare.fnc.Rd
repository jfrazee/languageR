\name{simulateLatinsquare.fnc}
\alias{simulateLatinsquare.fnc}
\title{Simulate simple Latin Square data and compare models}
\description{
  This function creates a user-specified number of simulated split-plot 
  datasets, and compares mixed-effects models with the by-subject anova.
}
\usage{
simulateLatinsquare.fnc(dat, with = TRUE, trial = 0, nruns = 100, 
   nsub = NA, nitem = NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame with the structure of the data set \code{latinsquare}.}
  \item{with}{Logical, if TRUE, effect of SOA built into the data.}
  \item{trial}{A number which, if nonzero, gives the magnitude of a
     learning or a fatigue effect.}
  \item{nruns}{A number indicating the required number of simulation runs.}
  \item{nsub}{A number for the number of subjects.}
  \item{nitem}{A number for the number of items.}
	\item{...}{other parameters to be passed through to plotting functions.}
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
\examples{\dontrun{
  data(latinsquare)
	\dontrun{
	library(lme4)
  simulateLatinsquare.fnc(latinsquare, nruns=100)
	}
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
