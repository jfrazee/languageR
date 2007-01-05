\name{subjects.latinsquare.fnc}
\alias{subjects.latinsquare.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{By-subject analysis of simple Latin Square data sets}
\description{
  This function is called by \code{simulateLatinsquare.fnc} for by-subject
  analysis of simulated Latin Square datasets.  It is not intended for
  independent use.
}
\usage{
subjects.latinsquare.fnc(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame with variables RT or RTsim , SOA, Subject, Item,
    Group and List, as in the \code{latinsquare} data set.}
}
\value{
  A list with components
  \item{p}{The p-value of the by-subject anova.}
  \item{data}{The input dataset.}
  \item{model}{The fitted model.}
}
\author{R. H. Baayen}
\seealso{See also \code{\link{simulateLatinsquare.fnc}}.}
\examples{\dontrun{
data(latinsquare)
subjects.latinsquare.fnc(latinsquare)$p
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
