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
   ilabel = NA, fun = NA, pred = NA, control = NA, ranefs = NA, n = 100, 
   intr = NA,  lockYlim = TRUE, addlines = FALSE, 
   withList = FALSE, cexsize = 0.5, linecolor = 1, addToExistingPlot = FALSE, 
   verbose = TRUE, ...)
}
\arguments{
  \item{model}{ a \code{LMM} or \code{GLMM} model object of class \code{lmerMod} }
  \item{xlabel}{ label for X-axis (if other than the variable name in the
    original model formula)}
  \item{xlabs}{ character vector with labels for X-axes in multipanel plot (if 
  other than the variable names in the original model formula); if used, 
  \code{xlabel} should not be specified}
  \item{ylabel}{ label for Y-axis (if other than the variable name of 
    the dependent variable in the original model formula) }
  \item{ylimit}{ range for vertical axis; if not specified, this range will 
    be chosen such that all data points across all subplots, including HPD intervals, will be accommodated }
  \item{ilabel}{ label for the interaction shown in the lower right-hand margin of the plot, overriding the original variable name in the model formula}
  \item{fun}{ a function to be applied for transforming the dependent variable, 
    if \code{NA}, no transformation is applied; for models with \code{family = "binomial"}, 
      fun is set to \code{plogis} by default; this can be disabled by setting
      \code{fun=function(x)return(x)}.}
  \item{pred}{ character string with name of predictor; 
    if specified, a single plot will produced for the partial effect of this
    specific predictor}
  \item{control}{ a two-element list \code{list(predictor, val)} specifying a \code{predictor} the value of which has to be set to \code{val} in the partial effect plot(s); the predictor name should be exactly as specified in \code{names(model@fixef)}.  It is up to the user to make sure that name and value make sense, the code here hands full 'control' to the user.} 
  \item{ranefs}{a four-element list \code{Group, Level, Predictor}, specifying a random-effect Group (e.g. \code{Subject}), a level (e.g., \code{S10}) and a value (e.g., \code{LogFrequency}) for which partial effects have to be calibrated.}
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
  \item{lockYlim}{ logical specifying whether all subplots should have the same
  range of values for the vertical axis; if \code{TRUE}, this range will be
  chosen to accomodate all fitted values including HDP intervals for all
  predictors across all plots}
  \item{addlines}{ if TRUE, adds line(s) between levels of same factor(s)}
  \item{withList}{ logical, if \code{TRUE}, a list will be output with all data 
  frames for the subplots}
  \item{cexsize}{ character expansion size (cex) for additional information in
  the plot for interactions}
  \item{linecolor}{ color of lines in the plot, by default set to 1 (black) }
  \item{addToExistingPlot}{ default FALSE, if set to TRUE, plot will be added to previous plot, but only if pred is specified}
  \item{verbose}{ if TRUE (default), effect sizes and default transformations are reported}
  \item{...}{ further graphical parameters to be passed down; warning: \code{col}, 
  \code{pch}, \code{lty} and \code{cex} will often generate an error as they are 
  internally already fully specified for specialized subplots}
}
\details{

  When no predictor is specified, a series of plots is produced for the partial
  effects of each predictor.  The graphs are shown for the reference level for
  factors and are adjusted for the median value for the other numerical
  predicors in the model.  Interactions are not shown.  The user should set up
  the appropriate number of subplots on the graphics device before running
  plotLMER.fnc().

  Instead of showing all predictors jointly, plotLMER.fnc() can also be used to
  plot the partial effect of a specific predictor.  When a specific predictor
  is specified (with \code{pred = ...}), a single plot is produced for that
  predictor.  In this case, the \code{intr} argument can be used to specify a
  single second predictor that enters into an interaction with the selected
  main predictor.  
  
  Polynomials have to be fitted with \code{poly(..., degree, raw=TRUE)} and
  restricted cubic splines with \code{rcs()} from the \code{rms} package.
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

###########################################################################
# we will stay as close to the older optimizer of lme4 as possible -
# this requires the optimx package and using the control option of lmer()
###########################################################################
require(optimx)

###########################################################################
# fitting a cosine with a spline (simulated data)
###########################################################################

require("rms", quietly=TRUE, character=TRUE)
require("lme4", quietly=TRUE, character=TRUE)
dfr = makeSplineData.fnc()
table(dfr$Subject)
xylowess.fnc(Y ~ X | Subject, data = dfr)
# the smoother doesn't recognize the cosine function implemented in makeSplineData.fnc()
dev.off()   

dfr.lmer = lmer(Y ~ rcs(X, 5) + (1|Subject), data = dfr,
  control = lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))
plotLMER.fnc(dfr.lmer)

# comparison with ols from Design package
dfr.lm = lm(Y~Subject+rcs(X), data=dfr, x=T, y=T)
dfr$fittedOLS = fitted(dfr.lm)
dfr$fittedLMER = as.vector(dfr.lmer@pp$X \%*\% fixef(dfr.lmer))

# we plot the lmer() fit in blue, the ols() fit in red (both adjusted for
# subject S1), and plot the underlying model in green

plot(dfr[dfr$Subject=="S1",]$X, 
  dfr[dfr$Subject=="S1",]$fittedLMER + ranef(dfr.lmer)[[1]]["S1",], 
  col="blue", ylim = c(24,30), xlab="X", ylab="Y", type="n")   

lines(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$fittedOLS, col="red")
lines(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$fittedLMER, col="blue")
lines(dfr[dfr$Subject=="S1",]$X, dfr[dfr$Subject=="S1",]$y+
  ranef(dfr.lmer)[[1]]["S1",], col="green")
legend(2,30,c("30+cos(x)", "lmer (S1)", "ols (S1)"), lty=rep(1,3), 
  col=c("green", "blue", "red"))


#############################################################
# a model with a raw polynomial
#############################################################

bg.lmer = lmer(LogRT ~ PC1+PC2+PC3 + ReadingScore +
  poly(OrthLength, 2, raw=TRUE) + LogFrequency + LogFamilySize +
  (1|Word) + (1|Subject)+(0+OrthLength|Subject) +
  (0+LogFrequency|Subject), data = beginningReaders,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

pars = par()
par(mfrow=c(3,3), mar=c(5,5,1,1))
plotLMER.fnc(bg.lmer, fun=exp, ylabel = "RT (ms)")

#############################################################
# a model with an interaction involving numeric predictors
#############################################################

danish.lmer = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
  ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
  poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
  (1|Subject) + (1|Word) + (1|Affix), data = danish,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))
danish.lmerA = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
  ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
  poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
  (1|Subject) + (1|Word) + (1|Affix), data = danish,
  subset=abs(scale(resid(danish.lmer)))<2.5,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

# plot for reference level of Sex
plotLMER.fnc(danish.lmerA, pred = "LogAffixFreq", 
  intr=list("LogWordFreq", round(quantile(danish$LogWordFreq),3), "beg",
  list(c("red", "green", "blue", "yellow", "purple"), rep(1,5))), 
  ylimit=c(6.5,7.0))

# this model has a significant three-way interaction
# for visualization, we can either relevel Sex and refit,
# or make use of the control option. First releveling:

danish$Sex=relevel(danish$Sex, "F")
danish.lmerF = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
  ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
  poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
  (1|Subject) + (1|Word) + (1|Affix), data = danish,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))
danish$Sex=relevel(danish$Sex, "M")
danish.lmerM = lmer(LogRT ~ PC1 + PC2 + PrevError + Rank +
  ResidSemRating + ResidFamSize + LogWordFreq*LogAffixFreq*Sex +  
  poly(LogCUP, 2, raw=TRUE) + LogUP + LogCUPtoEnd + 
  (1|Subject) + (1|Word) + (1|Affix), data = danish,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

# Next preparing for using the control option:
#
# names(fixef(danish.lmer))[10]     # SexM
# unique(danish.lmer@pp$X[,10])     # 1 0

par(mfrow=c(2,2))

plotLMER.fnc(danish.lmer, pred="LogWordFreq", ylimit=c(6.5,7.0),
intr=list("LogAffixFreq", round(quantile(danish$LogAffixFreq),2), "end"),
control=list("SexM", 0))
mtext("females", line=1.5, cex=0.9)

plotLMER.fnc(danish.lmer, pred="LogWordFreq", ylimit=c(6.5,7.0),
intr=list("LogAffixFreq", round(quantile(danish$LogAffixFreq),2), "end"),
control=list("SexM", 1))
mtext("males", line=1.5, cex=0.9)

plotLMER.fnc(danish.lmerF, pred="LogWordFreq", ylimit=c(6.5,7.0), 
intr=list("LogAffixFreq", round(quantile(danish$LogAffixFreq),2), "end"))
mtext("females", line=1.5, cex=0.9)

plotLMER.fnc(danish.lmerM, pred="LogWordFreq", ylimit=c(6.5, 7.0),
intr=list("LogAffixFreq", round(quantile(danish$LogAffixFreq),2), "end"))
mtext("males", line=1.5, cex=0.9)

par(mfrow=c(1,1))

#############################################################
# calculating effect sizes, defined as max - min
#############################################################

# effect size for a covariate

dfr = plotLMER.fnc(danish.lmerA, pred = "LogCUP", withList=TRUE)
max(dfr$LogCUP$Y)-min(dfr$LogCUP$Y)

# effect size for a factor

dfr = plotLMER.fnc(danish.lmerA, pred = "PrevError", withList=TRUE)
max(dfr$PrevError$Y)-min(dfr$PrevError$Y)


# effect sizes for the quantiles in an interaction plot

dfr = plotLMER.fnc(danish.lmerA, pred = "LogAffixFreq", 
  withList=TRUE,
  intr=list("LogWordFreq", round(quantile(danish$LogWordFreq),3), "beg"))

unlist(lapply(dfr$LogAffixFreq, FUN=function(X)return(max(X$Y)-min(X$Y))))


#############################################################
# plotting an interaction between two factors
#############################################################

danish$WordFreqFac = danish$LogWordFreq > median(danish$LogWordFreq)
danish.lmer2 = lmer(LogRT ~ WordFreqFac*Sex +  
  (1|Subject) + (1|Word) + (1|Affix), data = danish,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

plotLMER.fnc(danish.lmer2, pred = "Sex", 
  intr=list("WordFreqFac", c("TRUE", "FALSE"), "end", 
  list(c("red",  "blue"), rep(1,2))),
  ylimit=c(6.7,6.9), cexsize=1.0, addlines=TRUE)

#############################################################
# a generalized linear mixed-effects model
#############################################################

dative.lmer = glmer(RealizationOfRecipient ~ 
  AccessOfTheme + AccessOfRec + LengthOfRecipient + AnimacyOfRec +
  AnimacyOfTheme + PronomOfTheme + DefinOfTheme + LengthOfTheme +
  SemanticClass + Modality + (1|Verb), 
  data = dative, family = "binomial",
  control=glmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

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
