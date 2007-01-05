\name{items.split.plot.fnc}
\alias{items.split.plot.fnc}
\title{By-item anova for simulated split plot design}
\description{
  By-item anova for simulated data set as created within 
  simulateSplitPlot.fnc.   Not intended for independent use.
  Depends on the packages MASS, coda and lme4.
}
\usage{
items.split.plot.fnc(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Simulated data set with Subjects, Item, and SOA treatment,
  as created within simulateSplitPlot.fnc, or the splitplot dataset.}
}
\value{
  A list with components
  \item{p}{p-value of F-test for SOA.}
  \item{data}{the input data.}
  \item{model}{the fitted model.}
}
\author{R. H. Baayen}
\seealso{See also \code{\link{simulateSplitPlot.fnc}}.}
\examples{\dontrun{
data(splitplot)
items.split.plot.fnc(splitplot)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
