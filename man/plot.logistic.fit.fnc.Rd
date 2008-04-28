\name{plot.logistic.fit.fnc}
\alias{plot.logistic.fit.fnc}
\title{Plot for goodness of fit of logistic regression}
\description{
  This function plots observed proportions against mean predicted
  probabilities. For a good fit, points should be approximately on
  a straight line.
}
\usage{
plot.logistic.fit.fnc(x, data, method, where, ...)
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
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  A plot is produced on the graphics device.
}
\author{R. H. Baayen}
\examples{\dontrun{
data(dative)
library(lme4, keep.source=FALSE)

dative.glmm = lmer(RealizationOfRecipient ~ AccessOfTheme + 
   AccessOfRec + LengthOfRecipient + AnimacyOfRec +
   AnimacyOfTheme + PronomOfTheme + DefinOfTheme + LengthOfTheme +
   SemanticClass + Modality + (1|Verb), 
   data = dative, family = "binomial")

plot.logistic.fit.fnc (dative.glmm, dative)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
