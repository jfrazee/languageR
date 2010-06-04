\name{make.reg.fnc}
\alias{make.reg.fnc}
\title{Make a simulated data set with regression design}
\description{
  This convenience function creates a regression data set with subjects, items,
  and three numerical predictors, and optionally an effect of learning or
  fatigue.  This function is called by simulateRegression.fnc, and is not
  intended for independent use.
}
\usage{
make.reg.fnc(nsubj = 10, nitem = 20, beta = c(400, 2, 6, 4), 
learn = FALSE, learnRate = 10, stdevItem = 40, stdevSubj = 80, 
  stdevError = 50)
}
\arguments{
  \item{nsubj}{Number of subjects (random effect) required.}
  \item{nitem}{Number of items (random effect) required.}
  \item{beta}{A numeric vector with four beta weights: one for the 
     intercept and one for each of three predictors.}
  \item{learn}{A logical variable, if TRUE, a learning or fatigue effect
     will be implemented, as specified by \code{learnRate}.}
  \item{learnRate}{A number indicating learning (if negative) or fatigue
     (if positive).}
  \item{stdevItem}{A number specifying the standard deviation of the Item random effect.}
  \item{stdevSubj}{A number specifying the standard deviation of the Subject random effect.} 
  \item{stdevError}{A number specifying the standard deviation of the Residual Error.}
}
\value{
  A data frame with intercept, predictors labelled X, Y and Z, Item, Subject, the simulated
  random effects for Item and Subject, the residual errors, and the simulated RTs.
}
\author{R. H. Baayen}
\seealso{\code{\link{simulateRegression.fnc}}}
\examples{\dontrun{
  simdat = make.reg.fnc()
  library(lme4, keep.source = FALSE)
  simdat.lmer = lmer(RT ~ X + Y + Z + (1|Subject) + (1|Item), 
    data = simdat)
  simdat.lmer

  simdat = make.reg.fnc(learn = TRUE)
  simdat.lmer = lmer(RT ~ X + Y + Z + Trial + (1|Subject) + (1|Item), 
    data = simdat)
  simdat.lmer
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
