\name{spanish}
\alias{spanish}
\docType{data}
\title{Relative frequencies of tag trigrams is selected Spanish texts}
\description{
  Relative frequencies of the 120 most frequent tag trigrams in 15 texts
  contributed by 3 authors.
}
\usage{data(spanish)}
\format{
  A data frame with 120 observations on 15 variables documented in
  \code{spanishMeta}.
}
\references{
Spassova, M. S. (2006) \emph{Las marcas sintacticas de atribucion forense de
autoria de textos escritos en espanol}, Masters thesis, Institut Universitari
de Linguistica Aplicada, Universitat Pompeu Fabra, Barcelona. 
}
\examples{\dontrun{
data(spanish)
data(spanishMeta)

# principal components analysis

spanish.t = t(spanish)
spanish.pca = prcomp(spanish.t, center = TRUE, scale = TRUE)
spanish.x = data.frame(spanish.pca$x)
spanish.x = spanish.x[order(rownames(spanish.x)), ]

library(lattice)
splom(~spanish.x[ , 1:3], groups = spanishMeta$Author)

# linear discriminant analysis

library(MASS)
spanish.pca.lda = lda(spanish.x[ , 1:8], spanishMeta$Author)
plot(spanish.pca.lda)

# cross-validation

n = 8
spanish.t = spanish.t[order(rownames(spanish.t)), ]
predictedClasses = rep("", 15)
for (i in 1:15) {
  training = spanish.t[-i,]                     
  trainingAuthor = spanishMeta[-i,]$Author
  training.pca = prcomp(training, center=TRUE, scale=TRUE)
  training.x = data.frame(training.pca$x)
  training.x = training.x[order(rownames(training.x)), ]
  training.pca.lda = lda(training[ , 1:n], trainingAuthor)
  predictedClasses[i] = 
  as.character(predict(training.pca.lda, spanish.t[ , 1:n])$class[i])  
}

ncorrect = sum(predictedClasses==as.character(spanishMeta$Author))
ncorrect
sum(dbinom(ncorrect:15, 15, 1/3))
}}
\keyword{datasets}
