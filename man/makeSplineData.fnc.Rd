\name{makeSplineData.fnc}
\alias{makeSplineData.fnc}
\title{generate simulated data set with nonlinear function}
\description{
 creates a data set with \code{Y ~ 30+cos(X)} for 10 subjects,
 to compare restricted cubic spline in \code{lmer} with the
 spline of \code{ols}.
}
\usage{
makeSplineData.fnc(intr=0)
}
\arguments{
  \item{intr}{integer denoting type of data set: with 0 a data set with
    simple spline is made, with 1 a data set with a parallel interaction,
    and with 2 a data set with a crossed interaction.}
}
\details{
  Requires \code{Design} package to be attached.
}
\value{
  A data frame with as values:
  \item{y}{\code{y = 30 + cos(X)}}
  \item{X}{ranges from 2 to 8.28}
  \item{Subject}{random-effects factor with 10 levels}
  \item{Ranef}{subjects-specific changes to intercept}
  \item{Error}{by-observation noise}
  \item{Y}{the dependent variable, \code{y+Ranef+Error}}
}
\author{R. H. Baayen}
\note{ 
  intended for illustration only
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
\dontrun{
require("Design", quietly=TRUE, character=TRUE)
dfr = makeSplineData.fnc()
table(dfr$Subject)
xylowess.fnc(Y ~ X | Subject, data = dfr)

dfr.lmer = lmer(Y ~ rcs(X, 5) + (1|Subject), data = dfr)
dfr$fittedLMER = as.vector(dfr.lmer@X \%*\% fixef(dfr.lmer))

dfr.dd = datadist(dfr)
options(datadist='dfr.dd')
dfr.ols = ols(Y~Subject+rcs(X), data=dfr, x=T, y=T)
dfr$fittedOLS = fitted(dfr.ols)

# we plot the lmer() fit in blue, the ols() fit in red (both adjusted for
# subject S1), and plot the underlying model in green
plot(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$fittedLMER +
  ranef(dfr.lmer)[[1]]["S1",], type="l", col="blue",
  ylim = range(dfr$y + ranef(dfr.lmer)[[1]]["S1",],
  dfr[dfr$Subject == "S1",]$fittedLMER,
  dfr[dfr$Subject == "S1",]$fittedOLS), xlab="X", ylab="Y")   
lines(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$fittedOLS, col="red")
lines(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$y+ranef(dfr.lmer)[[1]]["S1",], 
  col="green")
legend(2,29,c("30+cos(x)", "lmer (S1)", "ols (S1)"), lty=rep(1,3), 
col=c("green", "blue", "red"))
}
}
\keyword{ regression }
