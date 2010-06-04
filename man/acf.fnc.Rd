\name{acf.fnc}
\alias{acf.fnc}
\title{Autocorrelation trellis graph}
\description{
  This function creates a trellis plot with autocorrelation functions
  for by-subject sequential dependencies in response latencies.
}
\usage{
acf.fnc(dat, group="Subject", time="Trial", x = "RT", plot=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame with (minimally) a grouping factor, an index for
   successive trails/events, and a behavioral measure}
  \item{group}{A grouping factor such as \code{Subject}}
  \item{time}{A sequential time measure such as \code{Trial} number in the
   experimental list}
  \item{x}{The dependent variable, usually a chronometric measure such as RT}
  \item{plot}{If true, a trellis graph is produced, otherwise a data frame
   with the data on which the trellis graph is based is returned}
  \item{...}{other optional arguments, such as \code{layout}}
}
\value{
  If \code{plot=TRUE}, a trellis graph, otherwise a data frame with as column 
  names
  \item{Lag}{Autocorrelation lag}
  \item{Acf}{Autocorrelation}
  \item{Subject}{The grouping factor, typically Subject}
  \item{ci}{The (approximate) 95\% confidence interval.}
}
\references{
R. H. Baayen (2001) \emph{Word Frequency Distributions}, Dordrecht: Kluwer.  
}
\author{R. H. Baayen}
%\seealso{}
\examples{\dontrun{
data(beginningReaders)
acf.fnc(beginningReaders, x="LogRT")   # autocorrelations even though nonword responses not included
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RT autocorrelations }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
