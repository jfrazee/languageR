\name{etymology}
\alias{etymology}
\docType{data}
\title{Etymological age and regularity in Dutch}
\description{
Estimated etymological age for regular and irregular monomorphemic Dutch
verbs, together with other distributional predictors of regularity.
}
\usage{data(etymology)}
\format{
  A data frame with 285 observations on the following 14 variables.
  \describe{
    \item{\code{Verb}}{a factor with the verbs as levels.}
    \item{\code{WrittenFrequency}}{a numeric vector of logarithmically transformed
         frequencies in written Dutch (as available in the CELEX lexical database).}
    \item{\code{NcountStem}}{a numeric vector for the number of orthographic neighbors.}
    \item{\code{MeanBigramFrequency}}{a numeric vector for mean log bigram frequency.}
    \item{\code{InflectionalEntropy}}{a numeric vector for Shannon's entropy calculated for the word's inflectional variants.}
    \item{\code{Auxiliary}}{a factor with levels \code{hebben}, \code{zijn} and \code{zijnheb} for the verb's auxiliary in the perfect tenses.}
    \item{\code{Regularity}}{a factor with levels \code{irregular} and \code{regular}.}
    \item{\code{LengthInLetters}}{a numeric vector of the word's orthographic length.}
    \item{\code{Denominative}}{a factor with levels \code{Den} and \code{N} specifying
       whether a verb is derived from a noun according to the CELEX lexical database.}
    \item{\code{FamilySize}}{a numeric vector for the number of types in the word's 
       morphological family.}
    \item{\code{EtymAge}}{an ordered factor with levels \code{Dutch}, \code{DutchGerman}, \code{WestGermanic}, \code{Germanic} and \code{IndoEuropean}.}
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
\dontrun{
data(etymology)

# ---- EtymAge should be an ordered factor, set contrasts accordingly

etymology$EtymAge = ordered(etymology$EtymAge, levels = c("Dutch",
"DutchGerman", "WestGermanic", "Germanic", "IndoEuropean")) 
options(contrasts=c("contr.treatment","contr.treatment"))

library(rms)
etymology.dd = datadist(etymology)
options(datadist = 'etymology.dd')

# ---- EtymAge as additional predictor for regularity

etymology.lrm = lrm(Regularity ~ WrittenFrequency + 
rcs(FamilySize, 3) + NcountStem + InflectionalEntropy + 
Auxiliary + Valency + NVratio + WrittenSpokenRatio + EtymAge, 
data = etymology, x = TRUE, y = TRUE)
anova(etymology.lrm)

# ---- EtymAge as dependent variable

etymology.lrm = lrm(EtymAge ~ WrittenFrequency + NcountStem +
MeanBigramFrequency + InflectionalEntropy + Auxiliary +
Regularity + LengthInLetters + Denominative + FamilySize + Valency + 
NVratio + WrittenSpokenRatio, data = etymology, x = TRUE, y = TRUE)

# ---- model simplification 

etymology.lrm = lrm(EtymAge ~ NcountStem + Regularity + Denominative, 
data = etymology, x = TRUE, y = TRUE)
validate(etymology.lrm, bw=TRUE, B=200)

# ---- plot partial effects and check assumptions ordinal regression

plot(Predict(etymology.lrm))
plot(etymology.lrm)
resid(etymology.lrm, 'score.binary', pl = TRUE)
plot.xmean.ordinaly(EtymAge ~ NcountStem, data = etymology)
}
}
\keyword{datasets}
