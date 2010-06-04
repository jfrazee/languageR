\name{oldFrench}
\alias{oldFrench}
\docType{data}
\title{Frequencies of tag trigrams in Old French texts}
\description{
Frequencies of 35 morphosyntactic tag trigrams in 343 Old French texts.
}
\usage{data(oldFrench)}
\format{
  A data frame with the frequencies of 35 tag trigrams (columns) for 343 Old French texts (rows)
  in the Nouveau Corpus d'Amsterdam.  See oldFrenchMeta for details on the texts (and manuscript
  versions).
}
\source{
  Data from Nouveau Corpus d'Amsterdam, \url{http://www.uni-stuttgart.de/lingrom/stein/corpus/}.
}
\references{
Ernestus, M., van Mulken, M. and Baayen, R. H. (2007)
De syntax van Oud-Franse ridders en heiligen in ruimte en tijd
To appear in \emph{Onze Taal}.
}
\examples{
data(oldFrench)
data(oldFrenchMeta)

oldFrench.ca = corres.fnc(oldFrench)

plot(oldFrench.ca, rlabels = oldFrenchMeta$Genre, 
rcol = as.numeric(oldFrenchMeta$Genre), rcex = 0.5, 
extreme = 0.1, ccol = "blue")
}
\keyword{datasets}
