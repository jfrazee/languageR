\name{plotLMER.fnc}
\alias{plotLMER.fnc}
\title{plot a mer object}
\description{
  Plot partial effects of a (generalized) linear mixed-effects model fit with
  \code{lmer}.  For gaussian models, 95\% highest posterior density credible 
  intervals can be added.
}
\usage{
plotLMER.fnc(model, xlabel = NA, xlabs = NA, ylabel = NA, ylimit = NA, 
   fun = NA, pred = NA, n = 100, intr = NA, mcmcMat = NA, lockYlim = TRUE, 
   addlines = FALSE, withList = FALSE, cexsize = 0.5, ...)
}
\arguments{
  \item{model}{ a \code{mer} model object }
  \item{xlabel}{ label for X-axis (if other than the variable name in the
    original model formula)}
  \item{xlabs}{ character vector with labels for X-axes in multipanel plot (if 
  other than the variable names in the original model formula); if used, 
  \code{xlabel} should not be specified}
  \item{ylabel}{ label for Y-axis (if other than the variable name of 
    the dependent variable in the original model formula) }
  \item{ylimit}{ range for vertical axis; if not specified, this range will 
    be chosen such that all data points across all subplots, including HPD intervals, will be accommodated }
  \item{fun}{ a function to be applied for transforming the dependent variable, 
    if \code{NA}, no transformation is applied; for models with \code{family = "binomial"}, 
      fun is set to \code{plogis} by default; this can be disabled by setting
      \code{fun=function(x)return(x)}.}
  \item{pred}{ character string with name of predictor; 
    if specified, a single plot will produced for the partial effect of this
    specific predictor}
  \item{n}{ integer denoting number of points for the plot, chosen at equally
    spaced intervals across the empirical range of the predictor variable }
  \item{intr}{ a list specifying an interaction to be graphed; obligatory
  arguments are (1) the name of the interaction variable, followed by (2) 
  a vector of values for that variable, followed by (3) the position for 
  interaction labels ('"beg"', '"mid"', or '"end"', or 'NA' if no labels are 
  desired), optionally followed by (4) a list with as first element
  a vector of colors and as second element a vector of line types. The number
  of elements in both vectors should match the number of values specified 
  under (2) for the interaction predictor.}
  \item{mcmcMat}{ an MCMC matrix as produced by \code{mcmcsamp} or 
  \code{pvals.fnc()$mcmc} for the input model; when this matrix is specified, 
  an attempt is made to add MCMC-based HPD intervals to the predicted values
  in the plot.  Works only for R 2.6.2 or higher, and for lme4 0.999 or higher. }
  \item{lockYlim}{ logical specifying whether all subplots should have the same
  range of values for the vertical axis; if \code{TRUE}, this range will be
  chosen to accomodate all fitted values including HDP intervals for all
  predictors across all plots}
  \item{addlines}{ if TRUE, adds line(s) between levels of same factor(s)}
  \item{withList}{ logical, if \code{TRUE}, a list will be output with all data 
  frames for the subplots}
  \item{cexsize}{ character expansion size (cex) for additional information in
  the plot for interactions}
  \item{...}{ further graphical parameters to be passed down; warning: \code{col}, 
  \code{pch}, \code{lty} and \code{cex} will often generate an error as they are 
  internally already fully specified for specialized subplots}
}
\details{

  When no predictor is specified, a series of plots is produced for the partial
  effects of each predictor.  The graphs are shown for the reference level for
  factors and are adjusted for the median value for the other numerical
  predicors in the model.  Interactions are not shown.  The user should set up the
  appropriate number of subplots on the graphics device before running plotLMER.fnc().

  Instead of showing all predictors jointly, plotLMER.fnc() can also be used to
  plot the partial effect of a specific predictor.  When a specific predictor
  is specified (with \code{pred = ...}), a single plot is produced for that
  predictor.  In this case, the \code{intr} argument can be used to specify a
  single second predictor that enters into an interaction with the selected
  main predictor.  
  
  Polynomials have to be fitted with \code{poly(..., degree, raw=TRUE)} and
  restricted cubic splines with \code{rcs()} from the \code{Design} package.

  When an MCMC object (as produced by the MCMC code in lme4, version 1.0 or
  higher) is specified (\code{mcmcMat = ... }), an attempt is made to plot
  MCMC HPD intervals for the predicted values.  This code is under development.
  For each MCMC sample, the posterior values of the fixed effects are used
  to predict the dependent variable for a range of its values, with other 
  predictors in the model held constant (factors at the reference level, 
  covariates at their medians).  For each of the predicted 'posterior' values,
  the HPD interval is obtained, and added to the graph.

}
\value{
  A plot is produced on the graphical device.
}
\references{The 'danish' dataset in the example section is contributed by Laura Winther-Balling, 
see Winther-Balling, L. and Baayen, R. H., Morphological effects in auditory word recognition: 
Evidence from Danish, Language and Cognitive Processes, in press.}
\author{R. H. Baayen}
\note{ 
  This code needs much more work, including (i) extension to \code{poly} with \code{raw=FALSE},
  and (ii) general clean-up of the code.
}
\seealso{See also other utilities in languageR for facilitating work with \code{lmer}}
\examples{
  \dontrun{

#############################################################
# fitting a cosine with a spline (simulated data)
#############################################################

require("Design", quietly=TRUE, character=TRUE)
dfr = makeSplineData.fnc()
table(dfr$Subject)
xylowess.fnc(Y ~ X | Subject, data = dfr)
# the smoother doesn't recognize the cosine function implemented in makeSplineData.fnc()
dev.off()   

dfr.lmer = lmer(Y ~ rcs(X, 5) + (1|Subject), data = dfr)
plotLMER.fnc(dfr.lmer)

# comparison with ols from Design package
dfr.dd = datadist(dfr)
options(datadist='dfr.dd')
dfr.ols = ols(Y~Subject+rcs(X), data=dfr, x=T, y=T)
dfr$fittedOLS = fitted(dfr.ols)
dfr$fittedLMER = as.vector(dfr.lmer@X \%*\% fixef(dfr.lmer))

# we plot the lmer() fit in blue, the ols() fit in red (both adjusted for
# subject S1), and plot the underlying model in green

plot(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$fittedLMER +
  ranef(dfr.lmer)[[1]]["S1",], type="l", col="blue",
  ylim = range(dfr$y + ranef(dfr.lmer)[[1]]["S1",],
  dfr[dfr$Subject == "S1",]$fittedLMER,
  dfr[dfr$Subject == "S1",]$fittedOLS),
  xlab="X", ylab="Y")   
lines(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$fittedOLS, 
col="red")
lines(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$y+
ranef(dfr.lmer)[[1]]["S1",], col="green")
legend(2,29,c("30+cos(x)", "lmer (S1)", "ols (S1)"), lty=rep(1,3), 
col=c("green", "blue", "red"))





#############################################################
# a model with a raw polynomial: the mcmc credible intervals
# will be different each time this example is run; for large
# numbers of samples, the plotting code will become quite slow
#############################################################

bg.lmer = lmer(LogRT ~ PC1+PC2+PC3 + ReadingScore +
poly(OrthLength, 2, raw=TRUE) + LogFrequency + LogFamilySize +
(1|Word) + (1|Subject)+(0+OrthLength|Subject) +
(0+LogFrequency|Subject), data = beginningReaders)
mcmc = pvals.fnc(bg.lmer, nsim=1000, withMCMC=TRUE)

pars = par()
par(mfrow=c(3,3), mar=c(5,5,1,1))
plotLMER.fnc(bg.lmer, mcmcMat=mcmc$mcmc, fun=exp, ylabel = "RT (ms)")

#############################################################
# a model with an interaction involving numeric predictors
#############################################################

danish.lmer = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
(1|Subject) + (1|Word) + (1|Affix), data = danish)
danish.lmerA = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
(1|Subject) + (1|Word) + (1|Affix), data = danish,
subset=abs(scale(resid(danish.lmer)))<2.5)
mcmc = pvals.fnc(danish.lmerA, nsim=10000, withMCMC=TRUE)
mcmc$fixed[,1:5]

# this model has a significant three-way interaction
# for visualization, we split the data by the sex of the subject

males = danish[abs(scale(resid(danish.lmer)))<2.5 & danish$Sex=="M",]
females = danish[abs(scale(resid(danish.lmer)))<2.5 & danish$Sex=="F",]

# we now fit individual models to these subsets of the data

males.lmer = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq +  
poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
(1|Subject) + (1|Word) + (1|Affix), data = males)
females.lmer = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq +  
poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
(1|Subject) + (1|Word) + (1|Affix), data = females)
femMCMC = pvals.fnc(females.lmer, nsim=500, withMCMC=TRUE, addPlot=FALSE)

# for which we then plot the interaction 
# of word frequency by affix frequency 

par(mfrow=c(2,2),mar=c(5,5,2,2))
plotLMER.fnc(females.lmer, pred = "LogAffixFreq", 
intr=list("LogWordFreq", round(quantile(females$LogWordFreq),3), "beg",
list(c("red", "green", "blue", "yellow", "purple"), rep(1,5))), ylimit=c(6.5,7.0),
mcmcMat=femMCMC$mcmc, xlabel = "log affix frequency", ylabel = "log RT auditory lexical decision")
mtext("females", line=0.5, cex=0.8)
plotLMER.fnc(males.lmer, pred = "LogAffixFreq", 
intr=list("LogWordFreq", round(quantile(males$LogWordFreq),3), "beg",
list(c("red", "green", "blue", "yellow", "purple"), 1:5)),
ylimit=c(6.5,7.0))
mtext("males", line=0.5, cex=0.8)
plotLMER.fnc(females.lmer, pred = "LogWordFreq", 
intr=list("LogAffixFreq", round(quantile(females$LogAffixFreq),3), "beg"),
ylimit=c(6.5,7.0), cexsize=1.2)
mtext("females", line=0.5, cex=0.8)
plotLMER.fnc(males.lmer, pred = "LogWordFreq", 
intr=list("LogAffixFreq", round(quantile(males$LogAffixFreq),3), "end"),
ylimit=c(6.5,7.0), cexsize=1.0)
mtext("males", line=0.5, cex=0.8)
par(mfrow=c(1,1))


#############################################################
# plotting an interaction between two factors
#############################################################

danish$WordFreqFac = danish$LogWordFreq > median(danish$LogWordFreq)
danish.lmer2 = lmer(LogRT ~ WordFreqFac*Sex +  
(1|Subject) + (1|Word) + (1|Affix), data = danish)
mcmc = pvals.fnc(danish.lmer2,nsim=1000,withMCMC=TRUE)

plotLMER.fnc(danish.lmer2, pred = "Sex", 
intr=list("WordFreqFac", c("TRUE", "FALSE"), "end", 
list(c("red",  "blue"), rep(1,2))),
ylimit=c(6.5,7.0), cexsize=1.0, addlines=TRUE, mcmcMat=mcmc$mcmc)

#############################################################
# a generalized linear mixed-effects model
#############################################################

dative.lmer = lmer(RealizationOfRecipient ~ 
AccessOfTheme + AccessOfRec + LengthOfRecipient + AnimacyOfRec +
AnimacyOfTheme + PronomOfTheme + DefinOfTheme + LengthOfTheme +
SemanticClass + Modality + (1|Verb), 
data = dative, family = "binomial")

par(mfrow=c(3,4),mar=c(5,5,1,1))
plotLMER.fnc(dative.lmer, fun=plogis, addlines=TRUE)

# with user-specified labels for the x-axis
par(mfrow=c(3,4),mar=c(5,5,1,1))
plotLMER.fnc(dative.lmer, fun=plogis, addlines=TRUE,
xlabs=unlist(strsplit("abcdefghij","")))

par(pars)


  }
}
\keyword{ regression }
