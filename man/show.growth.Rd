\name{show.growth}
\alias{show.growth}
\title{Plot method for growth objects.}
\description{
  A print method for growth objects created with \code{growth.fnc}.
}
\usage{
show.growth(x)
}
\arguments{
  \item{x}{A growth object.}
}
\value{
  Prints growth object.  To access the data frame embedded in the
  growth object, use \code{<my.growth.object>@data$data}.
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{growth.fnc}}.}
\examples{
  data(alice)
  alice.growth = growth.fnc(alice, chunks= c(5000, 10000, 15000))
  alice.growth
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
