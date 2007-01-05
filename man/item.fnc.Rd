\name{item.fnc}
\alias{item.fnc}
\title{Function for by-item regression used by simulateRegression.fnc}
\description{
  This function carries out a by-item regression for the simulated
  data sets generated in simulate.regression.fnc.  It is not designed
  to be used independently.
}
\usage{
item.fnc(data)
}
\arguments{
  \item{data}{A data frame as produced by make.reg.fnc().}
}
\value{
  A model fitted with lm().
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{simulateRegression.fnc}} and
\code{\link{make.reg.fnc}}.}
\examples{\dontrun{
  dat = make.reg.fnc()
  dat.lm = item.fnc(dat)
  summary(dat.lm)
}}
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
