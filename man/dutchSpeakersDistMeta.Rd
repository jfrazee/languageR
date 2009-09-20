\name{dutchSpeakersDistMeta}
\alias{dutchSpeakersDistMeta}
\docType{data}
\title{Metadata for dutchSpeakersDist}
\description{
Meta-data for the cross-entropy based between-speaker 
distance matrix \code{dutchSpeakersDist}
}
\usage{data(dutchSpeakersDistMeta)}
\format{
  A data frame with 165 observations on the following 6 variables.
  \describe{
    \item{\code{Speaker}}{a factor with speakers as levels.}
    \item{\code{Sex}}{a factor with levels \code{female} and \code{male}.}
    \item{\code{AgeYear}}{a numeric vector with the speakers' year of birth.}
    \item{\code{AgeGroup}}{a factor with levels \code{age18to24},
          \code{age25to34}, \code{age35to44}, \code{age45to55}, and
          \code{age56up}}. 
    \item{\code{ConversationType}}{a factor with levels \code{femaleOnly} \code{maleFemale}, \code{maleOnly}, and \code{unknown}.}
    \item{\code{EduLevel}}{a factor with levels \code{EduUnknown}, \code{high},
          \code{low} \code{mid}}
	}
}
\source{
  http://lands.let.kun.nl/cgn/
}
\references{
Juola, P. (2003) The time course of language change,
\emph{Computers and the Humanities}, 37, 77-96.

Juola, P. and Baayen, R. H. (2005) A Controlled-corpus Experiment in Authorship
Identification by Cross-entropy, \emph{Literary and Linguistic Computing}, 20,
59-67.
}
\examples{
	\dontrun{
     data(dutchSpeakersDistMeta)
  }
}
\keyword{datasets}
