\name{splitplot}
\alias{splitplot}
\docType{data}
\title{Simulated data set with split plot design}
\description{
  Simulated lexical decision latencies with priming as treatment and reaction
  time in lexical decision as dependent variable.
}
\usage{data(splitplot)}
\format{
  A data frame with 800 observations on the following 11 variables.
  \describe{
    \item{\code{items}}{A factor with levels \code{w1}, \code{w2},
      ...,  \code{w40}, coding 40 word items.}
    \item{\code{ritems}}{The by-word random adjustments to the intercept.}
    \item{\code{list}}{A factor with levels \code{listA} and \code{listB}. 
        The priming effect is counterbalanced for subjects across these 
        two lists, compare \code{table(splitplot$list, splitplot$subjects)}.}
    \item{\code{rlist}}{The by-list random adjustments to the intercept.}
    \item{\code{priming}}{A treatment factor with levels \code{primed} and
          \code{unprimed}.}
    \item{\code{fpriming}}{The priming effect, -30 for the primed and 0 for
          the unprimed condition.}
    \item{\code{subjects}}{A factor with levels \code{s1}, \code{s2}, 
          ... \code{s20} coding 20 subjects.}
    \item{\code{rsubject}}{The by-subject random adjustments to the intercept.}
    \item{\code{error}}{The by-observation noise.}
    \item{\code{int}}{The intercept.}
    \item{\code{RT}}{The reaction time.}
  }
}
\source{
R. H. Baayen, D. J. Davidson and D. M. Bates.  Mixed-effects modeling
with crossed random effects for subjects and items.  Manuscript under revision for \emph{Journal of Memory and Language}.
}
\examples{\dontrun{
data(splitplot)
table(splitplot$list, splitplot$subjects)
library(lme4)
dat.lmer1 = lmer2(RT ~ list*priming+(1+priming|subjects)+(1+list|items),data=dat)
dat.lmer2 = lmer2(RT ~ list*priming+(1+priming|subjects)+(1|items),data=dat)
dat.lmer3 = lmer2(RT ~ list*priming+(1|subjects)+(1|items),data=dat)
dat.lmer4 = lmer2(RT ~ list*priming+(1|subjects),data=dat)
anova(dat.lmer1, dat.lmer2, dat.lmer3, dat.lmer4)
print(dat.lmer3, corr=FALSE)
}}
\keyword{datasets}
