\name{summary.corres}
\alias{summary.corres}
\title{Summarize a correspondence object}
\description{
  This function provides a concise summary of a correspondence object.
}
\usage{
summary.corres(object, n = 2, returnList = FALSE, head = TRUE, ...)
}
\arguments{
  \item{object}{A correspondence object as produced by \code{corres}.}
  \item{n}{A number indicating number of dimensions to be summarized. }
  \item{returnList}{Logical, if TRUE, a list is returned with as components
    the full information on each factor, instead of only the first 6 lines.}
  \item{head}{Logical, if TRUE, first 6 rows of factor summaries are shown.}
	\item{...}{Additional arguments passed on to summaries.}
}
\value{
  A summary with eigenvalue rates, and coordinates, correlations,
  and contributions for the factors (by default, 2, unless n is set to 
  a higher number).  
}
\author{R. H. Baayen}
\seealso{See also \code{\link{corres.fnc}}.}
\examples{
  data(oldFrench)
  oldFrench.ca = corres.fnc(oldFrench)
  oldFrench.ca
  summary(oldFrench.ca)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
