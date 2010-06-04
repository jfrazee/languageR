\name{plot.logistic.fit.fnc}
\alias{plot.logistic.fit.fnc}
\title{Plot for goodness of fit of logistic regression}
\description{
  This function plots observed proportions against mean predicted
  probabilities. For a good fit, points should be approximately on
  a straight line.
}
\usage{
plot.logistic.fit.fnc(x, data, method, where, scalesize, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A logistic regression model fitted with \code{lmer} or
  \code{lrm}.}
  \item{data}{The data set to which the model was fitted.}
  \item{method}{Either "cut", in which case the vector of cut-off points
  supplied by the "where" argument will be used to partition the fitted 
  probabilities, or "shingle", in which a shingle (using \code{equal.count} 
  and its defaults) will be used.}
  \item{where}{A vector of cut-off points for partitioning the vector of
  fitted probabilities, by default \code{seq(0, 1, by=0.1)}}.
  \item{scalesize}{A positive real <= 1.  If not NA (the default), the circles
  representing data points in the graph are scaled to reflect the number of
  data points in the underlying data set.  The scalesize parameter specifies
  how large the largest circle will be compared to 1 inch.  For counts with
  large outliers, small values of scalesize are better.  See example below.
  }
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  A plot is produced on the graphics device.  The R-squared value shown
  above the plot represents the correlation between the X and Y values in the
  plot.  It does NOT represent the R-squared of the lrm or lmer model.
}
\author{R. H. Baayen}
\examples{\dontrun{
data(dative)
library(lme4, keep.source=FALSE)

dative.lrm = lrm(RealizationOfRecipient ~ AccessOfTheme + 
   AccessOfRec + LengthOfRecipient + AnimacyOfRec +
   AnimacyOfTheme + PronomOfTheme + DefinOfTheme + LengthOfTheme +
   SemanticClass + Modality, data = dative)

dative.glmm = lmer(RealizationOfRecipient ~ AccessOfTheme + 
   AccessOfRec + LengthOfRecipient + AnimacyOfRec +
   AnimacyOfTheme + PronomOfTheme + DefinOfTheme + LengthOfTheme +
   SemanticClass + Modality + (1|Verb), 
   data = dative, family = "binomial")

par(mfrow=c(2,2))
plot.logistic.fit.fnc (dative.lrm, dative)
mtext("lrm", 3, 3)
plot.logistic.fit.fnc (dative.glmm, dative)
mtext("lmer", 3, 3)
plot.logistic.fit.fnc (dative.lrm, dative, scalesize=0.2)
mtext("lrm", 3, 3)
plot.logistic.fit.fnc (dative.glmm, dative, method="shingle")
mtext("lmer", 3, 3)
par(mfrow=c(1,1))


}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
