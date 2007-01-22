\name{imaging}
\alias{imaging}
\docType{data}
\title{fMRI Filtered Signal and Priming Scores for Brain-Damaged Patients}
\description{
Filtered fMRI signal at the most significant voxel and average priming
scores for brain-damaged patients, in a study addressing the extent 
to which phonological and semantic processes recruit the same brain areas.
}
\usage{data(imaging)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{\code{Condition}}{a factor with levels 
       \code{irregulars} (the morphological condition involving priming using
       inflected forms of irregular English verbs, e.g., 'began'-'begin')
       and \code{semantics} (priming with semantically related words such 
       as 'card' and 'paper').}
    \item{\code{BehavioralScore}}{a numeric vector for the average priming
       scores.}
    \item{\code{FilteredSignal}}{a numeric vector for the intensity of the
       filtered fMRI signal at the most significant voxel.}
  }
}
\details{
  Location of data points reconstructed from the pixel map of Figure 2b
  of Tyler et al. 2005.
}
\source{
  Tyler, L.K.,  Marslen-Wilson, W.D. and Stamatakis, E.A. (2005)
  Differentiating lexical form, meaning, and structure
  in the neural language system, \emph{PNAS}, 102, 8375-8380.
}
\examples{
data(imaging)

imaging.lm = lm(FilteredSignal~BehavioralScore*Condition, data=imaging)
summary(imaging.lm)

plot(imaging$BehavioralScore, imaging$FilteredSignal, type = "n", 
  xlim = c(-30, 40), ylim = c(0, 80))
semantics = imaging[imaging$Condition == "semantics",]
irregulars = imaging[imaging$Condition == "irregulars",]
points(semantics$BehavioralScore, semantics$FilteredSignal, col = "black")
points(irregulars$BehavioralScore, irregulars$FilteredSignal, col = "darkgrey")
abline(lm(FilteredSignal ~ BehavioralScore, data = semantics), col = 'black')
abline(lm(FilteredSignal ~ BehavioralScore, data = irregulars), 
  col = 'darkgrey')

# model criticism

plot(imaging.lm)
outliers = c(1, 19) # given Cook's distance, or perhaps only
outliers = 1        # the outlier in the semantics subset
imaging.lm = lm(FilteredSignal ~ BehavioralScore * Condition, 
  data = imaging[-outliers, ])
summary(imaging.lm)


}
\keyword{datasets}
