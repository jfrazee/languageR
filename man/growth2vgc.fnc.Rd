\name{growth2vgc.fnc}
\alias{growth2vgc.fnc}
\title{Conversion of growth object into a vgc object}
\description{
This function converts a growth object (as defined in the languageR package)
to a vgc object (as defined in the zipfR package).
}
\usage{
growth2vgc.fnc(growth)
}
\arguments{
  \item{growth}{A growth object obtained with growth.fnc().}
}
\value{
  A vcg object as defined in the zipfR library.
}
\references{
  R. H. Baayen (2001) \emph{Word Frequency Distributions},
  Dordrecht: Kluwer Academic Publishers.

  zipfR Website: <URL: http://purl.org/stefan.evert/zipfR/>
}
\author{R. H. Baayen}
\seealso{See also \code{\link{growth.fnc}}and the zipfR package.}
\examples{\dontrun{
library(zipfR)

data(alice)
alice.growth = growth.fnc(text = alice, size = 648, nchunks = 40)
alice.vgc = growth2vgc.fnc(alice.growth)
plot(alice.vgc)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
