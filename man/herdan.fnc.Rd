\name{herdan.fnc}
\alias{herdan.fnc}
\title{Herdan's C}
\description{
  This function calculates Herdan's constant \code{C}.
}
\usage{
herdan.fnc(text, chunks)
}
\arguments{
  \item{text}{A vector of strings representing a text.}
  \item{chunks}{A vector of chunk sizes for which Herdan's C is required.
  Duplicate chunk sizes are not allowed, and the number of chunks should
  be at least 2.}
}
\value{
  A list with components
  \item{growth}{A data frame with token and type counts for the requested
     chunk sizes.}
  \item{C}{Herdan's C.}
}
\references{
Herdan, G. (1960) \emph{Type-Token Mathematics}, The Hague: Mouton.

Herdan, G. (1964) \emph{Quantitative Linguistics}, London: Buttersworths.
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{growth.fnc}}.}
\examples{
  data(alice)
  herdan.fnc(alice, cumsum(rep(floor(length(alice)/40), 40)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
