\name{aovlmer.fnc}
\alias{aovlmer.fnc}
\title{Compute p-values for factors in mixed models}
\description{
	This function computes p-values for factors in a mixed-effects model.  A
	p-value is returned based on an MCMC sample, as well as the anova table
	output by lmer extended with the p-values based on denominator degrees of
	freedom equal to the number of observations minus the number of fixed-effects
	coefficients.  For small datasets, this p-value is anticonservative.
}
\usage{
aovlmer.fnc(object, mcmc, which, noMCMC, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An lmer or glmer model for a response variable 
     fitted with \code{lmer}.}
  \item{mcmc}{A Markov chain Monte Carlo sample obtained with \code{mcmcsamp} 
	   for the lmer model.}
  \item{which}{A vector of integers or strings denoting the rows in the
	   table of coefficients that specify the coefficients for the factor.}
	\item{noMCMC}{A logical indicating whether MCMC sampling should be skipped,
	   in which case a data frame with a standard ANOVA table is returned.}
  \item{\dots}{Other optional arguments.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  When \code{noMCMC == FALSE}, a list with components
  \item{MCMC}{A list with the empirical p-value for the hypothesis 
	  that the columns of an MCMC sample specified by \code{which} have 
		mean zero versus a general multivariate distribution with elliptical 
		contours, and the rownames specified by \code{which}. 
	}
  \item{Ftests }{An anova table listing the (anticonservative) p-value 
	based on the F-test.  This table is also returned when noMCMC is set to TRUE.}
}
\author{R. H. Baayen, D. Bates}
\seealso{See also lmer, mcmcsamp.}
\examples{
	\dontrun{
	library(lme4, keep.source=FALSE)
	library(coda)
	data(latinsquare)
	l.lmer = lmer(RT~SOA+(1|Word)+(1|Subject), data=latinsquare)
	mcmc = mcmcsamp(l.lmer, n=10000)
	aovlmer.fnc(l.lmer, mcmc, c("SOAmedium", "SOAshort"))
	}
}
\keyword{ regression }
