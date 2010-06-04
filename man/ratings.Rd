\name{ratings}
\alias{ratings}
\docType{data}
\title{Ratings for 81 English nouns}
\description{
 Subjective frequency ratings, ratings of estimated weight, and ratings of estimated
 size, averaged over subjects, for 81 concrete English nouns.
}
\usage{data(ratings)}
\format{
  A data frame with 81 observations on the following 14 variables.
  \describe{
    \item{\code{Word}}{a factor with words as levels.} 
    \item{\code{Frequency}}{a numeric vector of logarithmically transformed frequencies}
    \item{\code{FamilySize}}{a numeric vector of logarithmically transformed 
                morphological family sizes.}
    \item{\code{SynsetCount}}{a numeric vector with logarithmically transformed counts
                of the number of synonym sets in WordNet in which the word is listed.}
    \item{\code{Length}}{a numeric vector for the length of the word in letters.}
    \item{\code{Class}}{a factor with levels \code{animal} and \code{plant}.}
    \item{\code{FreqSingular}}{a numeric vector for the frequency of the word in the
                singular.}
    \item{\code{FreqPlural}}{a numeric vector with the frequency of the word in the
                plural.}
    \item{\code{DerivEntropy}}{a numeric vector with the derivational entropies of the
                words.}
    \item{\code{Complex}}{a factor coding morphological complexity with levels 
               \code{complex} and \code{simplex}.}
    \item{\code{rInfl}}{a numeric vector coding the log of ratio of singular to plural
               frequencies.}
    \item{\code{meanWeightRating}}{a numeric vector for the estimated weight of the
               word's referent, averaged over subjects.}
    \item{\code{meanSizeRating}}{a numeric vector for the estimated size of the word's
               referent, averaged over subjects.}
    \item{\code{meanFamiliarity}}{a numeric vector with subjective frequency estimates,
               averaged over subjects.}
  }
}
\source{
Data collected together with Jen Hay at the University of Canterbury, Christchurch,
New Zealand, 2004.
}
\examples{\dontrun{
data(ratings)

ratings.lm = lm(meanSizeRating ~ meanFamiliarity * Class + 
I(meanFamiliarity^2), data = ratings)

ratings$fitted = fitted(ratings.lm)

plot(ratings$meanFamiliarity, ratings$meanSizeRating,       
xlab = "mean familiarity", ylab = "mean size rating", type = "n")
text(ratings$meanFamiliarity, ratings$meanSizeRating, 
substr(as.character(ratings$Class), 1, 1), col = 'darkgrey')

plants = ratings[ratings$Class == "plant", ]    
animals = ratings[ratings$Class == "animal", ]  
plants = plants[order(plants$meanFamiliarity),]
animals = animals[order(animals$meanFamiliarity),]

lines(plants$meanFamiliarity, plants$fitted)
lines(animals$meanFamiliarity, animals$fitted)
}}
\keyword{datasets}
