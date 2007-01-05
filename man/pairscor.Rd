\name{pairscor.fnc}
\alias{pairscor.fnc}
\title{Scatterplot matrix with correlations}
\description{
A matrix of scatterplots is produced with Pearson and Spearman correlations
in the lower triangle.  By default, smoothers are added to panels in the 
upper triangle, and histograms are added to the panels on the diagonal.  
}
\usage{
pairscor.fnc(data, hist = TRUE, smooth = TRUE,
  cex.points = 1,  col.points = "darkgrey")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or matrix with numeric vectors.}
  \item{hist}{a logical indicating whether panels on the diagonal
      should contain a histogram.}
  \item{smooth}{a logical indicating whether panels in the upper
      triangle should have a smoother added.}
  \item{cex.points}{a number indicating the size of the points in the
      panels in the upper triangle, available only when smoothers are added.}
  \item{col.points}{a number or string indicating the color of the points 
      in the panels in the upper triangle, available only when smoothers
      are added.}
}
\author{R. Harald Baayen}
\seealso{See Also \code{\link{pairs}} and \code{\link{panel.smooth}}.}
\examples{\dontrun{
  data(lexicalMeasures)
  pairscor.fnc(lexicalMeasures[,c("CelS", "Vf", "Ient", "NsyS", "Ncou")])
}}
%function(data) {
%	pairs(data, lower.panel = pairscor.lower, upper.panel = pairscor.upper)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
