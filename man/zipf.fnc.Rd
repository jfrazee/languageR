\name{zipf.fnc}
\alias{zipf.fnc}
\title{Zipf's rank frequency distribution}
\description{
  This function calculates Zipf's rank-frequency distribution for a text
  vector, and optionally produces the rank-frequency plot.
}
\usage{
zipf.fnc(text, plot = FALSE)
}
\arguments{
  \item{text}{A character vector containing a text.}
  \item{plot}{Logical, if TRUE, a rank-frequency plot is shown on the
              graphics device.}
}
\value{
  A data frame with variables
  \item{frequency }{Word Frequencies, ordered from large to small.}
  \item{freqOfFreq }{Frequencies of word frequencies.}
  \item{rank }{Zipf rank.}
}
\references{
Zipf, G. K. (1935) \emph{The Psycho-Biology of Language},
Boston: Houghton Mifflin.

Zipf, G. K. (1949) \emph{Human Behavior and the Principle of the Least Effort.
An Introduction to Human Ecology}, New York: Hafner.

Baayen, R. H. (2001) \emph{Word Frequency Distributions}, Dordrecht: Kluwer.
}
\author{R. H. Baayen}
\seealso{See also \code{\link{growth.fnc}}.}
\examples{
data(alice)
alice.zipf = zipf.fnc(alice, plot = TRUE) 
head(alice.zipf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
