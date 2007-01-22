\name{english}
\alias{english}
\docType{data}
\alias{english}
\title{English visual lexical decision and naming latencies}
\description{
  This data set gives mean visual lexical decision latencies and word
  naming latencies to 2284 monomorphemic English nouns and verbs,
  averaged for old and young subjects, with various predictor variables.
}
\usage{data(english)}
\format{
  A data frame with 4568 observations on the following variables.
  \describe{
    \item{\code{RTlexdec}}{numeric vector of log RT in visual lexical decision.}
    \item{\code{RTnaming}}{numeric vector of log RT in word naming.}
    \item{\code{Familiarity}}{numeric vector of subjective familiarity ratings.}
    \item{\code{Word}}{a factor with 2284 words.}
    \item{\code{AgeSubject}}{a factor with as levels the age group of the subject: 
          \code{young} versus \code{old}.} 
    \item{\code{WordCategory}}{a factor with as levels the word categories 
          \code{N} (noun) and \code{V} (verb).} 
    \item{\code{WrittenFrequency}}{numeric vector with log frequency in the CELEX lexical
          database.} 
    \item{\code{WrittenSpokenFrequencyRatio}}{numeric vector with the logged ratio 
          of written frequency (CELEX) to spoken frequency (British National Corpus).} 
    \item{\code{FamilySize}}{numeric vector with log morphological family size.} 
    \item{\code{DerivationalEntropy}}{numeric vector with derivational entropy.} 
    \item{\code{InflectionalEntropy}}{numeric vector with inflectional entropy.} 
    \item{\code{NumberSimplexSynsets}}{numeric vector with the log-transformed count of 
          synonym sets in WordNet in which the word is listed.} 
    \item{\code{NumberComplexSynsets}}{numeric vector with the log-transformed count of
          synonym sets in WordNet in which the word is listed as part of a compound.} 
    \item{\code{LengthInLetters}}{numeric vector with length of the word in letters.} 
    \item{\code{Ncount}}{numeric vector with orthographic neighborhood density.} 
    \item{\code{MeanBigramFrequency}}{numeric vector with mean log bigram frequency.} 
    \item{\code{FrequencyInitialDiphone}}{numeric vector with log frequency of initial diphone.} 
    \item{\code{ConspelV}}{numeric vector with type count of orthographic neighbors.} 
    \item{\code{ConspelN}}{numeric vector with token count of orthographic neighbors.} 
    \item{\code{ConphonV}}{numeric vector with type count of phonological neighbors.} 
    \item{\code{ConphonN}}{numeric vector with token count of phonological neighbors.} 
    \item{\code{ConfriendsV}}{numeric vector with type counts of consistent words.} 
    \item{\code{ConfriendsN}}{numeric vector with token counts of consistent words.} 
    \item{\code{ConffV}}{numeric vector with type count of forward inconsistent words }   
    \item{\code{ConffN}}{numeric vector with token count of forward inconsistent words  } 
    \item{\code{ConfbV}}{numeric vector with type count of backward inconsistent words  } 
    \item{\code{ConfbN}}{numeric vector with token count of backward inconsistent words } 
    \item{\code{NounFrequency}}{numeric vector with the frequency of the word used as noun.} 
    \item{\code{VerbFrequency}}{numeric vector with the frequency of the word used as verb.} 
    \item{\code{CV}}{factor specifying whether the initial phoneme of 
          the word is a consonant (\code{C}) or a vowel (\code{V}).} 
    \item{\code{Obstruent}}{factor specifying whether the initial phoneme 
          of the word is a continuant (\code{cont}) or an obstruent (\code{obst}).} 
    \item{\code{Frication}}{factor specifying whether the initial phoneme 
          has a burst (\code{burst}) or frication (\code{frication}) for 
          consonant-initial words, and for vowel-initial words whether the vowel is 
          \code{long} or \code{short}.} 
    \item{\code{Voice}}{factor indicating whether the initial phoneme is \code{voiced}
          or \code{voiceless}.} 
    \item{\code{FrequencyInitialDiphoneWord}}{numeric vector with the log-transformed 
          frequency of the initial diphone given that it is syllable-initial.} 
    \item{\code{FrequencyInitialDiphoneSyllable}}{numeric vector with the log-transformed 
          frequency of the initial diphone given that it is word initial.} 
    \item{\code{CorrectLexdec}}{numeric vector with the proportion of 
          subjects that accepted the item as a word in lexical decision.}
  }
}
\source{ 
Balota, D.A., Cortese, M.J. and Pilotti, M. (1999) \emph{
Visual lexical decision latencies for 2906 words}. Available at
\url{http://www.artsci.wustl.edu/~dbalota/lexical_decision.html}.

Spieler, D. H. and Balota, D. A. (1998) \emph{
Naming latencies for 2820 words},
\url{http://www.artsci.wustl.edu/~dbalota/naming.html}.
}
\references{ 
Balota, D., Cortese, M., Sergent-Marshall, S., Spieler, D. and Yap, M.
  (2004) Visual word recognition for single-syllable words, \emph{Journal of
  Experimental Psychology:General}, 133, 283-316.

Baayen, R.H., Feldman, L. and Schreuder, R. (2006)
   Morphological influences on the recognition of monosyllabic 
   monomorphemic words, \emph{Journal of Memory and Language},
   53, 496-512.
}
\examples{
data(english)

# ---- orthogonalize orthographic consistency measures

items = english[english$AgeSubject == "young",]
items.pca = prcomp(items[ , c(18:27)], center = TRUE, scale = TRUE)
x = as.data.frame(items.pca$rotation[,1:4])
items$PC1 =  items.pca$x[,1]
items$PC2 =  items.pca$x[,2]
items$PC3 =  items.pca$x[,3]
items$PC4 =  items.pca$x[,4]
items2 = english[english$AgeSubject != "young", ]
items2$PC1 =  items.pca$x[,1]
items2$PC2 =  items.pca$x[,2]
items2$PC3 =  items.pca$x[,3]
items2$PC4 =  items.pca$x[,4]
english = rbind(items, items2) 

# ---- add Noun-Verb frequency ratio

english$NVratio = log(english$NounFrequency+1)-log(english$VerbFrequency+1)

# ---- build model with ols() from Design

library(Design)
english.dd = datadist(english)
options(datadist = 'english.dd')

english.ols = ols(RTlexdec ~ Voice + PC1 + MeanBigramFrequency + 
   rcs(WrittenFrequency, 5) + rcs(WrittenSpokenFrequencyRatio, 3) + 
   NVratio + WordCategory + AgeSubject +
   rcs(FamilySize, 3) + InflectionalEntropy + 
   NumberComplexSynsets + rcs(WrittenFrequency, 5) : AgeSubject,
   data = english, x = TRUE, y = TRUE)

# ---- plot partial effects

par(mfrow = c(4, 3), mar = c(4, 4, 1, 1), oma = rep(1, 4))
plot(english.ols, adj.subtitle = FALSE, ylim = c(6.4, 6.9), conf.int = FALSE)
par(mfrow = c(1, 1))

# ---- validate the model

validate(english.ols, bw = TRUE, B = 200)

}
\keyword{datasets}
