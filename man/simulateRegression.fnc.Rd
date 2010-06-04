\name{simulateRegression.fnc}
\alias{simulateRegression.fnc}
\title{Simulate regression data and compare models}
\description{
  This function creates a user-specified number of simulated regression 
  datasets, and compares mixed-effects regression with random regression,
  by-subject regression, by-item regression, and by-subject plus by-item
  regression.  Optionally, an effect of learning or fatigue can be
  incorporated.  
}
\usage{
simulateRegression.fnc(beta = c(400, 2, 6, 4), nitem = 20, nsubj = 10, 
stdevItem = 40, stdevSubj = 80, stdevError = 50, nruns = 100, learn = FALSE, 
learnRate = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{A numeric vector with beta weights for the intercept and three
     numeric predictors.}
  \item{nitem}{A number specifying the number of items.}
  \item{nsubj}{A number specifying the number of subjects.}
  \item{stdevItem}{A number specifying the standard deviation of the
        Item random effect.}
  \item{stdevSubj}{ A number specifying the standard deviation of the
        Subject random effect.}
  \item{stdevError}{A number specifying the standard deviation of the
        Residual Error.}
  \item{nruns}{A number specifying the required number of simulated datasets.}
  \item{learn}{A logical that if TRUE, allows an effect of learning or
        fatigue into the model.}
  \item{learnRate}{A number specifying the learning rate (if negative) or
        the effect of fatigue (if positive).}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  A list with components
  \item{alpha05}{proportion of runs in which predictors are significant at the
        05 significance level.}
  \item{alpha01 }{proportion of runs in which predictors are significant at the
        01 significance level.}
  \item{ranef}{mean estimated random effects.}
  As this may take some time, the index of each completed run is shown
  on the output device.
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{make.reg.fnc}}.}
\examples{\dontrun{
	library(lme4)
  simulateRegression.fnc(beta = c(400, 2, 6, 4), nruns = 5)

  \dontrun{simulateRegression.fnc(beta = c(400, 2, 6, 0), nruns = 1000, learn = TRUE)}
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
