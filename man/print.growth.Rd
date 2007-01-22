\name{print.growth}
\alias{print.growth}
\title{Print method for growth objects.}
\description{
  Print method for growth objects.
}
\usage{
print.growth(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A growth object, as produced by \code{growth.fnc}.}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  The data frame with chunk sizes and associated vocabulary statistics
  is printed.   To access the data frame that is being shown, use
  \code{<my.growth.object>@data$data}.
}
\author{R. H. Baayen}
\seealso{See also \code{\link{growth.fnc}}.}
\examples{
  data(alice)
  alice.growth = growth.fnc(alice)
  alice.growth
  # for accessing the printed data frame:
  alice.growth@data$data[1:4,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
