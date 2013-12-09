\name{pvals.fnc}
\alias{pvals.fnc}
\title{Compute p-values and MCMC confidence intervals for mixed models}
\description{
  This function used to calculate p-values and HPD intervals for the 
	parameters of models fitted with \code{lmer}.  

  As MCMC is no longer supported by lme4, this function is now obsolete
  and does no longer produce any output, other than a warning.

  See the lme4 function pvalues() for alternatives.
}
\usage{
pvals.fnc(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{LMM} or \code{GLMM} model object of class \code{lmerMod}}
	\item{...}{Optional arguments that can be passed down.}
}
\value{
  A warning.
}
\author{R. H. Baayen}
\seealso{pvalues}
\examples{
  \dontrun{
  data(primingHeid) 
  library(lme4)

  # remove extreme outliers
  primingHeid = primingHeid[primingHeid$RT < 7.1,]

  # fit mixed-effects model

  # we will stay as close to the older optimizer of lme4 as possible -
  # this requires the optimx package and using the control option of lmer()

  require(optimx)
  require(lmerTest)

  primingHeid.lmer = lmer(RT ~ RTtoPrime * ResponseToPrime + 
    Condition + (1|Subject) + (1|Word), data = primingHeid,
    control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))
  summary(primingHeid.lmer)
  anova(primingHeid.lmer)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
