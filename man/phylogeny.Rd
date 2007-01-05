\name{phylogeny}
\alias{phylogeny}
\docType{data}
\title{Phylogenetic relations between Papuan and Oceanic languages}
\description{
Phylogenetic relations between Papuan and Oceanic languages: 127 
grammatical traits (absent/present) for 31 languages.
}
\usage{data(phylogeny)}
\format{
  A data frame with 31 observations on the following 127 variables.
  \describe{
    \item{\code{Language}}{a factor for 31 languages: \code{Anem}, \code{Ata},
       \code{Bali}, \code{Banoni}, \code{Bilua}, \code{Buin}, \code{Gapapaiwa},
       \code{Kairiru}, \code{Kaulong}, \code{Kilivila}, \code{Kokota}, 
       \code{Kol}, \code{Kuot}, \code{Lavukaleve}, \code{Mali}, \code{Motuna}, 
       \code{Nalik}, \code{Nasioi}, \code{Rotokas}, \code{Roviana}, 
       \code{Savosavo}, \code{Siar}, \code{Sisiqa}, \code{Sudest}, 
       \code{Sulka}, \code{Taiof}, \code{Takia}, \code{Touo}, \code{Tungag},
       \code{Yabem} and \code{Yeli_Dnye}.}
    \item{\code{Family}}{a factor with levels \code{Oceanic} and \code{Papuan}.}
    \item{\code{Frics}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PrenasalizedStops}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PhonDistBetweenLAndR}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PhonVelarFricOrGlide}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PhonVoicingContrAmongStops}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PhonConsLength}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PhonVowelLength}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ContrPhonTypesForVowels}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PhonStress}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{WordFinalConss}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ConsClusters}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{DefOrSpecArt}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{IndefOrNonSpecArt}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ArticleNounOrder}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NounInitNps}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{InclExclDist}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PronNum}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PronRelationship}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PronConflation}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{MoreThan2DegreesDistDem}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NonSpkrAnchoredDem}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerticalityDem}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ClassifiedDem}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NumDeterminedDecl}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{GenderDeterminedDecl}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SuppletiveNouns}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SingMarkedNoun}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{DualMarkedNoun}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PlMarkedNoun}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{OtherNumMarkedNoun}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{LimitedDistNumMarking}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NounClassesGenders}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ConcordBeyondNp}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NumeralClassifiers}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PossClassifiers}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PossClasses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{Inalienability}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{MultiplePossConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PrefixMarkedPoss}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SuffixMarkedPoss}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{MarkedPossr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{MarkedPossessee}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PossrPossdOrder}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{DecimalNumerals}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{QuinaryNumerals}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{CollectiveNouns}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{AdjVerbLexOverlap}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{AdjAttributionPred}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{CoreCaseMarking}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ObliqueCaseMarking}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{Prepositions}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{Postpositions}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{TamPerson}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbPrefixesProclitics}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbSuffixesEnclitics}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PunctualContinuous}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{RealisIrrealis}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SSuffix}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SPrefix}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ASuffix}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{APrefix}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{OSuffix}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{OPrefix}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbVarTam}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbVarVClass}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbVarClauseType}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbVarPerson}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NumStemAlt}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PersonStemAlt}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SepVerbNumPerson}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{Portmanteau3Plus}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{DistributedCategory}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NonCore}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{RecipientObj}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{X3PlacePreds}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbNeg}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbDirection}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbSuppletion}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ConjugationClasses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{TransIntransAlt}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{TransitivizingMorph}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{IntranstivizingMorph}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ReflexiveMorph}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ReciprocalMorph}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbClassifiers}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{Copula}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NonVbPreds}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SerialVerbConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{Auxiliaries}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbCompounds}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbAdjunctConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VbIncorporation}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ExistentialVerb}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{IrregularGive}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ClosedClassOfVb}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SvIntransClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VsIntransClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VInitTransClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VMedialTransClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VFinalTransClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{FixedConstituentOrder}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ClauseFinalNeg}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ClauseInitNeg}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ImpVs.DeclNeg}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VbAndNonVbPredIdentity}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SOMorphInBasicConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SAMorphInBasicConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SOMorphInComplexConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SAMorphInComplexConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SynConflationOfSO}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ControlledUncontrolled}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ClauseChaining}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SimultaneousSequential}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{SayInDesidConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{RelativeClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{PurpSubClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{TemporalSubClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{ComplementClauses}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{CausBySerialVerbConstr}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{CausByBoundAffClit}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{CausByConstrInvolvingSay}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{MorphTopicOrFocus}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{TailHeadLinkage}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{VerbRedup}}{a numeric vector, 1: presence, 0: absence}
    \item{\code{NounRedup}}{a numeric vector, 1: presence, 0: absence}
  }
}
\source{
Dunn, M., Terrill, A., Reesink, G., Foley, R. A. and Levinson, S. C. (2005)
Structural phylogenetics and the reconstruction of ancient language history,
\emph{Science}, 309, 2072-2075.
}
\examples{\dontrun{
data(phylogeny)
library(ape)

# joint analysis of both language families

phylogeny.dist = dist(phylogeny[ ,3:ncol(phylogeny)], method = "binary")
phylogeny.dist.tr = nj(phylogeny.dist)
families =
  as.character(phylogeny$Family[as.numeric(phylogeny.dist.tr$tip.label)])
languages =
  as.character(phylogeny$Language[as.numeric(phylogeny.dist.tr$tip.label)])
phylogeny.dist.tr$tip.label = languages
plot(phylogeny.dist.tr, type="u", font = as.numeric(as.factor(families)))

# analysis of subset of Papuan languages

papuan = phylogeny[phylogeny$Family == "Papuan",]
papuan$Language = as.factor(as.character(papuan$Language))
papuan.meta = papuan[ ,1:2]
papuan.mat = papuan[, 3:ncol(papuan)]
papuan.meta$Geography = c(
  "Bougainville", "Bismarck Archipelago", "Bougainville", 
  "Bismarck Archipelago", "Bismarck Archipelago", "Central Solomons",
  "Bougainville", "Louisiade Archipelago", "Bougainville", 
  "Bismarck Archipelago", "Bismarck Archipelago", 
  "Bismarck Archipelago", "Central Solomons", "Central Solomons", 
  "Central Solomons")
papuan.dist = dist(papuan.mat, method = "binary")
papuan.dist.tr = nj(papuan.dist)
fonts = 
  as.character(papuan.meta$Geography[as.numeric(
    papuan.dist.tr$tip.label)])
papuan.dist.tr$tip.label = 
  as.character(papuan.meta$Language[as.numeric(
    papuan.dist.tr$tip.label)])
plot(papuan.dist.tr, type = "u", font = as.numeric(as.factor(fonts)))
}}
\keyword{datasets}
