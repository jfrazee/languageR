\name{spanishFunctionWords}
\alias{spanishFunctionWords}
\docType{data}
\title{Relative frequencies of function words in selected Spanish texts}
\description{
  Relative frequencies of the 120 most frequent function words in 15 texts
  contributed by 3 authors.
}
\usage{data(spanishFunctionWords)}
\format{
   A data frame with 120 observations on 15 variables documented in
   \code{spanishMeta}.
}
\references{
Spassova, M. S. (2006) \emph{Las marcas sintacticas de atribucion forense de
autoria de textos escritos en espanol}, Masters thesis, Institut Universitari
de Linguistica Aplicada, Universitat Pompeu Fabra, Barcelona. 
}
\examples{
data(spanishFunctionWords)
data(spanishMeta)

# principal components analysis

spanishFunctionWords.t = t(spanishFunctionWords)
spanishFunctionWords.t = 
  spanishFunctionWords.t[order(rownames(spanishFunctionWords.t)), ]
spanishFunctionWords.pca = 
  prcomp(spanishFunctionWords.t, center = TRUE, scale = TRUE)

sdevs = spanishFunctionWords.pca$sdev^2
n = sum(sdevs/sum(sdevs)> 0.05) 

# linear discriminant analysis with cross-validation

library(MASS)

predictedClasses = rep("", 15)
for (i in 1:15) {
  training = spanishFunctionWords.t[-i,]                   
  trainingAuthor = spanishMeta[-i,]$Author
  training.pca = prcomp(training, center = TRUE, scale = TRUE)
  training.x = data.frame(training.pca$x)
  training.x = training.x[order(rownames(training.x)), ]
  training.pca.lda = lda(training[ , 1:n], trainingAuthor)
  cl = predict(training.pca.lda, spanishFunctionWords.t[,1:n])$class[i]
  predictedClasses[i] = as.character(cl)
}

ncorrect = sum(predictedClasses==spanishMeta$Author)
sum(dbinom(ncorrect:15, 15, 1/3))
}
\keyword{datasets}
