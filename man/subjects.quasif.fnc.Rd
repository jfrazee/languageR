\name{subjects.quasif.fnc}
\alias{subjects.quasif.fnc}
\title{By-subject analysis of data sets requiring quasi-F ratios}
\description{
  This function is called by \code{simulateQuasif.fnc} for by-subject
  analysis of simulated datasets traditionally requiring quasi-F ratios.  
	It is not intended for independent use.
}
\usage{
subjects.quasif.fnc(dat)
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
\seealso{See also \code{\link{simulateQuasif.fnc}}.}
\examples{\dontrun{
  data(quasif)
  subjects.quasif.fnc(quasif)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
