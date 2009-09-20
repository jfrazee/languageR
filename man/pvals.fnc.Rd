\name{pvals.fnc}
\alias{pvals.fnc}
\title{Compute p-values and MCMC confidence intervals for mixed models}
\description{
  This function calculates p-values and HPD intervals for the 
	parameters of models fitted with \code{lmer}.  For the fixed-effects
	parameters, anticonservative p-values based on the t statistic with
	the upper bound for the degrees of freedom are added.
}
\usage{
pvals.fnc(object, nsim = 10000, ndigits = 4, withMCMC = FALSE, 
addPlot=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An lmer or glmer model object 
     fitted with \code{lmer}.  Currently, MCMC sampling for generalized 
     linear mixed models may not work off the shelf when there is more 
     than one random intercept in the model.  The code for random effects
     structure with correlation parameters is not yet implemented in 
     \code{mcmcsamp}, this currently blocks calculation of the fixed and
     random subtables. }
  \item{nsim}{An integer denoting the required number of 
     Markov chain Monte Carlo samples. }
  \item{ndigits}{An integer denoting the number of 
     decimal digits in the output. }
  \item{withMCMC}{A logical indicating whether the output of \code{mcmcsamp}
     should be returned.}
  \item{addPlot}{If TRUE (default), a trellis graph for the posterior 
     distributions of the parameters is shown on the current
     graphics device.}
	\item{...}{Optional arguments that can be passed down.}
}
\value{
  A list with components
  \item{fixed}{A data frame with the estimated coefficients, their MCMC
	   mean, the HPD 95% confidence intervals, the MCMC probability,
		 and the probability based on the t distribution with the number of 
		 observations minus the number of fixed-effects coefficients as degrees 
		 of freedom.  This last p-value is anti-conservative, especially
		 for small data sets.
    }
  \item{random}{A table with the MCMC means and HPD 95% confidence intervals
	   for the standard deviations and correlations of the random effects.
     The HPD intervals are based on the untransformed values (log for standard
     deviations, atanhs for correlations), the endpoints of the HPD interval are
     back-transformed to the scale of the input data frame for ease of interpretation.
    }
  \item{mcmc}{The output of \code{mcmcsamp}.}
}
\author{R. H. Baayen}
\seealso{See also lmer, mcmcsamp, and HPDinterval.}
\examples{
  \dontrun{
  data(primingHeid) 
  library(lme4, keep.source=FALSE)

  # remove extreme outliers
  primingHeid = primingHeid[primingHeid$RT < 7.1,]

  # fit mixed-effects model

  primingHeid.lmer = lmer(RT ~ RTtoPrime * ResponseToPrime + 
    Condition + (1|Subject) + (1|Word), data = primingHeid)

  mcmc = pvals.fnc(primingHeid.lmer, nsim=10000, withMCMC=TRUE)

  mcmc$fixed
  mcmc$random
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
