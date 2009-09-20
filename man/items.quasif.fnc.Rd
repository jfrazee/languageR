\name{items.quasif.fnc}
\alias{items.quasif.fnc}
\title{By-item anova for simulated data for quasi-F analysis}
\description{
  By-item anova for simulated data set as created within 
  simulateQuasif.fnc.   Not intended for independent use.
  Depends on the packages MASS, coda and lme4.
}
\usage{
items.quasif.fnc(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Simulated data set with Subjects, Item, and SOA treatment,
  as created within simulateQuasif.fnc, or the quasif dataset.}
}
\value{
  A list with components
  \item{p}{p-value of F-test for SOA.}
  \item{data}{the input data.}
  \item{model}{the fitted model.}
}
\author{R. H. Baayen}
\seealso{See also \code{\link{simulateQuasif.fnc}}.}
\examples{\dontrun{
data(quasif)
items.quasif.fnc(quasif)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
