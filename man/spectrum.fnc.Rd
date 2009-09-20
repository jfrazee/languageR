\name{spectrum.fnc}
\alias{spectrum.fnc}
\title{Frequency spectrum from text vector}
\description{
  This function creates a frequency spectrum for a 
  text in character vector form.
}
\usage{
spectrum.fnc(text)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{text}{A character vector containing the words of a text.}
}
\value{
  A data frame with as column variables
  \item{frequency}{Word frequencies.}
  \item{freqOfFreq}{The frequencies of the word frequencies.}
}
\references{
R. H. Baayen (2001) \emph{Word Frequency Distributions}, Dordrecht: Kluwer.  
}
\author{R. H. Baayen}
\seealso{See Also the \code{zipfR} package.}
\examples{\dontrun{
data(alice)
alice.spectrum = spectrum.fnc(alice)
head(alice.spectrum)
tail(alice.spectrum)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
