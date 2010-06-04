\name{mvrnormplot.fnc}
\alias{mvrnormplot.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatterplot of bivariate standard normal distribution}
\description{
This function produces a scatterplot for a bivariate standard normal	
distribution with least squares regression line.
}
\usage{
mvrnormplot.fnc(r, n, limits)
}
\arguments{
  \item{r}{The correlation, defaults to 0.9.}
  \item{n}{Number of simulated data points, defaults to 100.}
  \item{limits}{Optional range for the axes}
}
\value{
  A scatterplot with ordinary least squares regression line is shown on the
	graphics device, with sample estimate of r added at the top of the plot.
}
\author{R. H. Baayen}
\seealso{mvrnorm (MASS package)}
\examples{\dontrun{
	mvrnormplot.fnc(r=0.9, n=100)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
