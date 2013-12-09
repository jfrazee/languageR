\name{shrinkage}
\alias{shrinkage}
\docType{data}
\title{Data set illustrating shrinkage}
\description{
  Simulated data set for illustrating shrinkage.
}
\usage{data(shrinkage)}
\format{
  A data frame with 200 observations on the following 6 variables.
  \describe{
    \item{\code{intercept}}{a numeric vector for the intercept.}
    \item{\code{frequency}}{a numeric vector for word frequency.}
    \item{\code{subject}}{a factor for subjects with levels \code{S1},
      \code{S2}, ... , \code{S10}.}
    \item{\code{error}}{a numeric vector for residuals.}
    \item{\code{ranef}}{a numeric vector for random effect.}
    \item{\code{RT}}{a numeric vector for simulated RTs.}
  }
}
\examples{\dontrun{
data(shrinkage)

require(lme4)
require(lmerTest)
require(optimx)

shrinkage.lmer = lmer(RT ~ frequency + (1|subject), 
  data = shrinkage,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb"))
shrinkage.lmList = lmList(RT ~ frequency | subject, data = shrinkage)

# and visualize the difference between random regression 
# and mixed-effects regression

mixed = coef(shrinkage.lmer)[[1]]
random = coef(shrinkage.lmList)
subj = unique(shrinkage[,c("subject", "ranef")])
subj = subj[order(subj$subject),]
subj$random = random[,1]
subj$mixed = mixed[,1]
subj = subj[order(subj$random),]
subj$rank = 1:nrow(subj)

par(mfrow=c(1,2))
plot(subj$rank, subj$random, xlab="rank", ylab="RT", ylim=c(200,550), type="n")
text(subj$rank, subj$random, as.character(subj$subject), cex=0.8, col="red")
mtext("random regression", 3, 1)
points(subj$rank, 400+subj$ranef, col="blue")
abline(h=400)
plot(subj$rank, subj$mixed, xlab="rank", ylab="RT", ylim=c(200,550), type="n")
text(subj$rank, subj$mixed, as.character(subj$subject), cex=0.8, col = "red")
mtext("mixed-effects regression", 3, 1)
points(subj$rank, 400+subj$ranef, col="blue")
abline(h=400)
par(mfrow=c(1,1))
}}
\keyword{datasets}
