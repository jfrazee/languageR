\name{print.corres}
\alias{print.corres}
\title{Print method for correspondence object}
\description{
  Prints eigenvalues and eigenvalue rates for a correspondence object.
}
\usage{
\method{print}{corres}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A correspondence object.}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  Report of eigenvalues and eigenvalue rates.
}
\author{R. H. Baayen}
\seealso{See also \code{\link{corres.fnc}}. }
\examples{\dontrun{
  data(oldFrench)
  oldFrench.ca = corres.fnc(oldFrench)
  oldFrench.ca
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
