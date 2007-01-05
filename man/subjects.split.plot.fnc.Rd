\name{subjects.split.plot.fnc}
\alias{subjects.split.plot.fnc}
\title{By-subject analysis of simple split-plot data sets}
\description{
  This function is called by \code{simulateSplitPlot.fnc} for by-subject
  analysis of simulated split-plot datasets.  It is not intended for
  independent use.
}
\usage{
subjects.split.plot.fnc(dat)
}
\arguments{
  \item{dat}{A data frame with variables RT or RTsim , SOA, Subject and Item.}
}
\value{
  A list with components
  \item{p}{p-value for by-subject F-test.}
  \item{data}{Data set with aggregated subject means.}
  \item{model}{Anova table of fitted model.}
}
\author{R. H. Baayen}
\seealso{See also \code{\link{simulateSplitPlot.fnc}}.}
\examples{\dontrun{
  data(splitplot)
  subjects.split.plot.fnc(splitplot)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
