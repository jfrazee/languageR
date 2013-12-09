\name{aovlmer.fnc}
\alias{aovlmer.fnc}
\title{Compute p-values for factors in mixed models}
\description{
  This function no longer works with recent versions of lme4.  
  For p-values, see the anova() function in the lmerTest package.
  }
\usage{
aovlmer.fnc(object,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An lmer or glmer model for a response variable 
     fitted with \code{lmer}.}
  \item{\dots}{Other optional arguments.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A warning message.
}
\author{R. H. Baayen, D. Bates}
\seealso{See anova in lmerTest.}
\examples{
  \dontrun{
  library(optimx)
  library(lme4)
  data(latinsquare)
  l.lmer = lmer(RT~SOA+(1|Word)+(1|Subject), data=latinsquare,
    control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))
  library(lmerTest)
  summary(l.lmer)
  anova(l.lmer)
  }
}
\keyword{ regression }
