\name{text2spc.fnc}
\alias{text2spc.fnc}
\title{Create a frequency spectrum from a text vector}
\description{
  This functions takes a text in the form of a character vector as input,
  and outputs a frequency spectrum object as defined in the \code{zipfR} package.
}
\usage{
text2spc.fnc(text)
}
\arguments{
  \item{text}{A text in the form of a character vector.}
}
\value{
  A \code{spc} spectrum object as defined in the \code{zipfR} package.
}
\author{R. H. Baayen}
\seealso{See the documentation for \code{zipfR} for \code{spc} objects.}
\examples{\dontrun{
  library(zipfR)
  data(alice)
  alice.spc = text2spc.fnc(alice)
  plot(alice.spc)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
