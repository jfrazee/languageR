\name{regularity}
\alias{regularity}
\docType{data}
\title{Regular and irregular Dutch verbs}
\description{
Regular and irregular Dutch verbs and selected lexical and distributional
properties.
}
\usage{data(regularity)}
\format{
  A data frame with 700 observations on the following 13 variables.
  \describe{
    \item{\code{Verb}}{a factor with the verbs as levels.}
    \item{\code{WrittenFrequency}}{a numeric vector of logarithmically transformed
         frequencies in written Dutch (as available in the CELEX lexical database).}
    \item{\code{NcountStem}}{a numeric vector for the number of orthographic neighbors.}
    \item{\code{VerbalSynsets}}{a numeric vector for the number of verbal synsets in WordNet.}
    \item{\code{MeanBigramFrequency}}{a numeric vector for mean log bigram frequency.}
    \item{\code{InflectionalEntropy}}{a numeric vector for Shannon's entropy calculated for the word's inflectional variants.}
    \item{\code{Auxiliary}}{a factor with levels \code{hebben}, \code{zijn} and \code{zijnheb} for the verb's auxiliary in the perfect tenses.}
    \item{\code{Regularity}}{a factor with levels \code{irregular} and \code{regular}.}
    \item{\code{LengthInLetters}}{a numeric vector of the word's orthographic length.}
    \item{\code{FamilySize}}{a numeric vector for the number of types in the word's 
       morphological family.}
    \item{\code{Valency}}{a numeric vector for the verb's valency, estimated by its 
        number of argument structures.}
    \item{\code{NVratio}}{a numeric vector for the log-transformed ratio of the nominal
        and verbal frequencies of use.}
    \item{\code{WrittenSpokenRatio}}{a numeric vector for the log-transformed ratio of the frequencies in written and spoken Dutch.}
  }
}
\references{
Baayen, R. H. and Moscoso del Prado Martin, F. (2005) Semantic density and
past-tense formation in three Germanic languages, Language, 81, 666-698.

Tabak, W., Schreuder, R. and Baayen, R. H. (2005) Lexical statistics and
lexical processing: semantic density, information complexity, sex, and
irregularity in Dutch, in Kepser, S. and Reis, M., \emph{Linguistic Evidence -
Empirical, Theoretical, and Computational Perspectives}, Berlin: Mouton de
Gruyter, pp. 529-555.
}
\examples{
data(regularity)

# ---- predicting regularity with a logistic regression model

library(Design)
regularity.dd = datadist(regularity)
options(datadist = 'regularity.dd')

regularity.lrm = lrm(Regularity ~ WrittenFrequency + 
rcs(FamilySize, 3) + NcountStem + InflectionalEntropy + 
Auxiliary + Valency + NVratio + WrittenSpokenRatio, 
data = regularity, x = TRUE, y = TRUE)

anova(regularity.lrm)

# ---- model validation

validate(regularity.lrm, bw = TRUE, B = 200)
pentrace(regularity.lrm, seq(0, 0.8, by = 0.05))
regularity.lrm.pen = update(regularity.lrm, penalty = 0.6)
regularity.lrm.pen

# ---- a plot of the partial effects

par(mfrow = c(3, 3))
plot(regularity.lrm.pen, fun = plogis, ylab = "Pr(regular)", 
adj.subtitle = FALSE, ylim = c(0, 1))
par(mfrow = c(1, 1))


# predicting regularity with a support vector machine

library(e1071)
regularity$AuxNum = as.numeric(regularity$Auxiliary)
regularity.svm = svm(regularity[, -c(1,8,10)], regularity$Regularity, cross=10)
summary(regularity.svm)
}
\keyword{datasets}
