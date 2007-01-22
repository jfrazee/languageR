\name{variationLijk}
\alias{variationLijk}
\docType{data}
\title{Variation in spoken Dutch in the use of the suffix -lijk}
\description{
  This dataset documents variation in the use of the suffix \emph{-lijk}, as
  realized in 32 words, in spoken Dutch across region (Flanders versus The
  Netherlands), sex (females versus males) and education (high versus mid).
}
\usage{data(variationLijk)}
\format{
  A data frame with 32 observations on the following 8 variables.
  \describe{
    \item{\code{nlfemaleHigh}}{a numeric vector with counts for Dutch 
       females with a mid education level.}
    \item{\code{nlfemaleMid}}{a numeric vector counts for Dutch 
       females with a high education level.}
    \item{\code{nlmaleHigh}}{a numeric vector counts for Dutch 
       males with a high education level.}
    \item{\code{nlmaleMid}}{a numeric vector counts for Dutch 
       males with a mid education level.}
    \item{\code{vlfemaleHigh}}{a numeric vector counts for Flemish 
       females with a high education level.}
    \item{\code{vlfemaleMid}}{a numeric vector counts for Flemish 
       females with a mid education level.}
    \item{\code{vlmaleHigh}}{a numeric vector counts for Flemish 
       males with a high education level.}
    \item{\code{vlmaleMid}}{a numeric vector counts for Flemish 
       males with a mid education level.}
  }
}
\references{
Keune, K., Ernestus, M., Van Hout, R. and Baayen, R.H. (2005) Social,
geographical, and register variation in Dutch: From written 'mogelijk' to
spoken 'mok', \emph{Corpus Linguistics and Linguistic Theory}, 1, 183-223.
}
\examples{
data(variationLijk)
variationLijk.ca = corres.fnc(variationLijk)
plot(variationLijk.ca, rcex=0.7, ccol="black",
  rcol = rep("blue", nrow(variationLijk))) 
}
\keyword{datasets}
