\name{beginningReaders}
\alias{beginningReaders}
\docType{data}
\title{Visual lexical decision with beginning readers}
\description{
Visual lexical decision latencies for beginning readers (8 year-old
Dutch children).
}
\usage{data(beginningReaders)}
\format{
  A data frame with 7923 observations on the following 13 variables.
  \describe{
    \item{\code{Word}}{a factor for the words.}
    \item{\code{Subject}}{a factor for the subjects.}
    \item{\code{LogRT}}{a numeric vector with the log-transformed reaction time
          (in ms).}
    \item{\code{Trial}}{a numeric vector coding the rank of the trial in the 
          experimental list.}
    \item{\code{OrthLength}}{a numeric vector coding the word's length in letters.}
    \item{\code{LogFrequency}}{a numeric vector with log-transformed frequency in
          Vermeer's frequency dictionary of Dutch children's texts.}
    \item{\code{LogFamilySize}}{a numeric vector with the log-transformed morphological
          family size count (with family members judged to be unknown to young children
          removed).}
    \item{\code{ReadingScore}}{a numeric vector with a score for reading proficiency.}
    \item{\code{ProportionOfErrors}}{a numeric vector for the proportion of error responses for the word.}
    \item{\code{PC1}}{a numeric vector for the first principal component of a PCA
          orthogonalization of the preceding 4 reaction times}
    \item{\code{PC2}}{a numeric vector for the second principal component of a PCA
          orthogonalization of the preceding 4 reaction times}
    \item{\code{PC3}}{a numeric vector for the third principal component of a PCA
          orthogonalization of the preceding 4 reaction times}
    \item{\code{PC4}}{a numeric vector for the fourth principal component of a PCA
          orthogonalization of the preceding 4 reaction times}
  }
}
\references{
  Perdijk, K., Schreuder, R., Verhoeven, L. and Baayen, R. H. (2006)
   \emph{Tracing individual differences in reading skills of young children with 
   linear mixed-effects models}.  Manuscript, Radboud University Nijmegen.
}
\examples{
\dontrun{
data(beginningReaders)
require(lme4)
require(optimx)
require(lmerTest)

beginningReaders.lmer = lmer(LogRT ~  PC1 + PC2 + PC3  + ReadingScore +
  OrthLength + I(OrthLength^2) + LogFrequency + LogFamilySize +
  (1|Word) + (1|Subject) + (0+LogFrequency|Subject) + 
  (0+OrthLength|Subject) + (0+PC1|Subject), 
  data = beginningReaders,
  control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))
summary(beginningReaders.lmer)
}
}
\keyword{datasets}
