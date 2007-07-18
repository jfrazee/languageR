\name{growth.fnc}
\alias{growth.fnc}
\title{Calculate vocabulary growth curve and vocabulary richness measures}
\description{
This function calculates, for an increasing sequence of text sizes,
the observed number of types, hapax legomena, dis legomena, tris legomena,
and selected measures of lexical richness.
}
\usage{
growth.fnc(text = alice, size = 646, nchunks = 40, chunks = 0)
}
\arguments{
  \item{text}{A vector of strings representing a text.}
  \item{size}{An integer giving the size of a text chunk 
     when the text is to be split into a series of equally-sized text chunks.}
  \item{nchunks}{An integer denoting the number of desired equally-sized
     text chunks.}
  \item{chunks}{An integer vector denoting the token sizes for which growth
     measures are required. When chunks is specified, \code{size} and 
     \code{nchunks} are ignored.}
}
\value{
  A growth object with methods for plotting, printing.
  As running this function on large texts may take some time,
  a period is printed on the output device for each completed chunk
  to indicate progress.

  The data frame with the actual measures, which can be extracted with 
  \code{object.name@data$data}, has the following columns.
  
    \item{\code{Chunk}}{a numeric vector with chunk numbers.} 
    \item{\code{Tokens}}{a numeric vector with the number of tokens up to
       and including the current chunk.} 
    \item{\code{Types}}{a numeric vector with the number of types up to and
       including the current chunk.} 
    \item{\code{HapaxLegomena}}{a numeric vector with the corresponding count
       of hapax legomena.} 
    \item{\code{DisLegomena}}{a numeric vector with the corresponding count
       of dis legomena.} 
    \item{\code{TrisLegomena}}{a numeric vector with the corresponding count
       of tris legomena.} 
    \item{\code{Yule}}{a numeric vector with Yule's \code{K}.}
    \item{\code{Zipf}}{a numeric vector with the slope of Zipf's rank-frequency
       curve in the double-logarithmic plane.}
    \item{\code{TypeTokenRatio}}{a numeric vector with the ratio of types to
       tokens.}
    \item{\code{Herdan}}{a numeric vector with Herdan's \code{C}.}
    \item{\code{Guiraud}}{a numeric vector with Guiraud's \code{R}.}
    \item{\code{Sichel}}{a numeric vector with Sichel's \code{S}.}
    \item{\code{Lognormal}}{a numeric vector with mean log frequency.}
}
\references{
   R. H. Baayen (2001) \emph{Word Frequency Distributions},
   Dordrecht: Kluwer Academic Publishers.

   Tweedie, F. J. & Baayen, R. H. (1998) How variable may a constant be?
   Measures of lexical richness in perspective, \emph{Computers and the
   Humanities}, 32, 323-352.
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{plot.growth}}, and the zipfR package.}
\examples{\dontrun{
  data(alice)
  alice.growth = growth.fnc(alice)
  plot(alice.growth)
}}
% R documentation directory.
\keyword{ models }
