\name{lexicalMeasures}
\alias{lexicalMeasures}
\docType{data}
\title{Lexical measures for 2233 English monomorphemic words}
\description{
  Lexical distributional measures for 2233 English monomorphemic words.  This dataset provides
  a subset of the data available in the dataset \code{english}.
}
\usage{data(lexicalMeasures)}
\format{
  A data frame with 2233 observations on the following 24 variables.
  \describe{
    \item{\code{Word}}{a factor with 2284 words.}
    \item{\code{CelS}}{numeric vector with log-transformed lemma frequency in the CELEX lexical
          database.} 
    \item{\code{Fdif}}{numeric vector with the logged ratio 
          of written frequency (CELEX) to spoken frequency (British National Corpus).} 
    \item{\code{Vf}}{numeric vector with log morphological family size.} 
    \item{\code{Dent}}{numeric vector with derivational entropy.} 
    \item{\code{Ient}}{numeric vector with inflectional entropy.} 
    \item{\code{NsyS}}{numeric vector with the log-transformed count of 
          synonym sets in WordNet in which the word is listed.} 
    \item{\code{NsyC}}{numeric vector with the log-transformed count of
          synonym sets in WordNet in which the word is listed as part of a compound.} 
    \item{\code{Len}}{numeric vector with length of the word in letters.} 
    \item{\code{Ncou}}{numeric vector with orthographic neighborhood density.} 
    \item{\code{Bigr}}{numeric vector with mean log bigram frequency.} 
    \item{\code{InBi}}{numeric vector with log frequency of initial diphone.} 
    \item{\code{spelV}}{numeric vector with type count of orthographic neighbors.} 
    \item{\code{spelN}}{numeric vector with token count of orthographic neighbors.} 
    \item{\code{phonV}}{numeric vector with type count of phonological neighbors.} 
    \item{\code{phonN}}{numeric vector with token count of phonological neighbors.} 
    \item{\code{friendsV}}{numeric vector with type counts of consistent words.} 
    \item{\code{friendsN}}{numeric vector with token counts of consistent words.} 
    \item{\code{ffV}}{numeric vector with type count of forward inconsistent words.}   
    \item{\code{ffN}}{numeric vector with token count of forward inconsistent words.} 
    \item{\code{fbV}}{numeric vector with type count of backward inconsistent words.} 
    \item{\code{fbN}}{numeric vector with token count of backward inconsistent words} 
    \item{\code{ffNonzero}}{a numeric vector with the count of forward inconsistent words 
       with nonzero frequency.}
    \item{\code{NVratio}}{a numeric vector with the logarithmically transformed ratio
        of the noun and verb frequencies.}
  }
}
\references{
Baayen, R.H., Feldman, L. and Schreuder, R. (2006)
   Morphological influences on the recognition of monosyllabic 
   monomorphemic words, \emph{Journal of Memory and Language},
   53, 496-512.
}
\examples{\dontrun{
data(lexicalMeasures)
data(lexicalMeasuresDist)

library(rms)
library(cluster)
plot(varclus(as.matrix(lexicalMeasures[,-1])))

lexicalMeasures.cor = cor(lexicalMeasures[,-1], method = "spearman")^2
lexicalMeasures.dist = dist(lexicalMeasures.cor)
pltree(diana(lexicalMeasures.dist))

data(lexicalMeasuresClasses)
x = data.frame(measure = rownames(lexicalMeasures.cor), 
cluster = cutree(diana(lexicalMeasures.dist), 5),
class = lexicalMeasuresClasses$Class)
x = x[order(x$cluster), ]
x
}}
\keyword{datasets}
