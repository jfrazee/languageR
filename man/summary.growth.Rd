\name{summary.growth}
\alias{summary.growth}
\title{Summary method for growth objects}
\description{
  Summary method for vocabulary growth objects created with \code{growth.fnc}.
}
\usage{
summary.growth(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A vocabulary growth object.}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  The growth object is printed. For access to the data frame inside
  the object, use \code{<my.growth.object>@data$data}.
}
\author{R. H. Baayen}
\seealso{See also \code{\link{growth.fnc}}.}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
