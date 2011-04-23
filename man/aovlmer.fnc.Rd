\name{aovlmer.fnc}
\alias{aovlmer.fnc}
\title{Compute p-values for factors in mixed models}
\description{
	This function computes p-values for factors in a mixed-effects model on
  the basis of an MCMC sample.  
  }
\usage{
aovlmer.fnc(object, mcmc, which, noMCMC, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An lmer or glmer model for a response variable 
     fitted with \code{lmer}.}
  \item{mcmc}{A Markov chain Monte Carlo sample obtained with 
     \code{pvals.fnc(..., withMCMC=TRUE)$mcmc} for the lmer model.}
  \item{which}{A vector of integers or strings denoting the rows in the
	   table of coefficients that specify the coefficients for the factor.}
	\item{noMCMC}{In earlier versions of \code{lmer}, a single-argument  \code{anova} was supported.  In more recent versions, this is no longer possible.  Hence the p-values based on F-tests are no longer supported.}
  \item{\dots}{Other optional arguments.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list with components
  \item{p}{ the empirical p-value for the hypothesis 
	  that the columns of an MCMC sample specified by \code{which} have 
		mean zero versus a general multivariate distribution with elliptical 
		contours}
  \item{which}{ the rownames specified by \code{which}}. 
}
\author{R. H. Baayen, D. Bates}
\seealso{See also lmer, mcmcsamp.}
\examples{
	\dontrun{
	library(lme4, keep.source=FALSE)
	data(latinsquare)
	l.lmer = lmer(RT~SOA+(1|Word)+(1|Subject), data=latinsquare)
  mcmc = pvals.fnc(l.lmer, nsim=10000, withMCMC=TRUE)
  aovlmer.fnc(l.lmer, mcmc$mcmc, c("SOAmedium", "SOAshort"))
	}
}
\keyword{ regression }
