\name{lexicalMeasuresClasses}
\alias{lexicalMeasuresClasses}
\docType{data}
\title{Classification of lexical measures}
\description{
  A data frame labelling the lexical measures in the 
  dataset \code{lexicalMeasures} as measures of form or meaning.
}
\usage{data(lexicalMeasuresClasses)}
\format{
  A data frame with 23 observations on the following 3 variables.
  \describe{
    \item{\code{Variable}}{a factor with as levels the measures:
      \describe{
        \item{\code{Bigr}}{Mean Bigram Frequency.}
        \item{\code{CelS}}{CELEX Frequency.} 
        \item{\code{Dent}}{Derivational Entropy.} 
        \item{\code{fbN}}{Token Count of Backward Inconsistent Words.} 
        \item{\code{fbV}}{Type Count of Backward Inconsistent Words.} 
        \item{\code{Fdif}}{Ratio of Frequencies in Written and Spoken English.}
        \item{\code{ffN}}{Token Count of Forward Inconsistent Words.} 
        \item{\code{ffNonzero}}{Type Count of Forward Inconsistent Words 
              with Nonzero Frequency.}
        \item{\code{ffV}}{Type Count of Forward Inconsistent Words}
        \item{\code{friendsN}}{Token Count of Consistent Words.}
        \item{\code{friendsV}}{Type Count of Consistent Words.}
        \item{\code{Ient}}{Inflectional Entropy}
        \item{\code{InBi}}{Initial Bigram Frequency}
        \item{\code{Len}}{Length in Letters}
        \item{\code{Ncou}}{Orthographic Neighborhood Density}
        \item{\code{NsyC}}{Number of Complex Synsets}
        \item{\code{NsyS}}{Number of Simplex Synsets}
        \item{\code{NVratio}}{Ratio of Noun and Verb Frequencies}
        \item{\code{phonN}}{Token Count of Phonological Neighbors.}
        \item{\code{phonV}}{Type Count of Phonological Neighbors.}
        \item{\code{spelN}}{Token Count of Orthographic Neighbors.}
        \item{\code{spelV}}{Type Count of Orthographic Neighbors.}
        \item{\code{Vf}}{Morphological Family Size.}
      }
    }
    \item{\code{Class}}{a factor with levels \code{Form} and \code{Meaning}.}
    \item{\code{Explanation}}{a factor with glosses for the variables.}
  }
}
\references{
Baayen, R.H., Feldman, L. and Schreuder, R. (2006)
   Morphological influences on the recognition of monosyllabic 
   monomorphemic words, \emph{Journal of Memory and Language},
   53, 496-512.
}
\examples{\dontrun{
library(cluster)
data(lexicalMeasures)
data(lexicalMeasuresClasses)

lexicalMeasures.cor = cor(lexicalMeasures[,-1], method = "spearman")^2
x = data.frame(measure = rownames(lexicalMeasures.cor), 
cluster = cutree(diana(dist(lexicalMeasures.cor)), 5),
class = lexicalMeasuresClasses$Class)
x = x[order(x$cluster), ]
x
}}
\keyword{datasets}
