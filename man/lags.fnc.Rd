\name{lags.fnc}
\alias{lags.fnc}
\title{Calculate vector at specified lag}
\description{
  This function calculates for a given dependent variable the value of
  that variable at lag timesteps earlier in the time series of an
  experiment.
}
\usage{
lags.fnc(dat, time="Trial", group = "Subject", depvar = "RT", lag=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame with (minimally) a grouping factor, an time index for
   successive trails/events, and a behavioral measure}
  \item{group}{A grouping factor such as \code{Subject}}
  \item{time}{A sequential time index measure such as \code{Trial} number in an
   experimental list}
  \item{depvar}{The dependent variable, usually a chronometric measure such as RT}
  \item{lag}{The lag for which previous values are to be extracted}
}
\value{
  A vector with the values of the dependent variable at the specified lag.
  The by-group mean is substituted for the first lag timestep(s), 
  for which there is/are no preceding value(s) for the dependent variable.
}
\author{R. H. Baayen}
%\seealso{}
\examples{\dontrun{
  dfr = data.frame(Subject=c(rep("a", 5), rep("b", 5)),
                   Trial = c(rep(1:5,2)),
                   RT = rnorm(10, 500, 40))
  dfr$prevRT = lag.fnc(dfr, time="Trial", group="Subject", depvar="RT")
  dfr

}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{autocorrelation lag}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
