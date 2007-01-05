\name{dutchSpeakersDist}
\alias{dutchSpeakersDist}
\docType{data}
\title{Cross-entropy based distances between speakers}
\description{
A distance matrix for the conversations of 165 speakers in the 
Spoken Dutch Corpus.  Metadata on the speakers are available in
a separate dataset, \code{dutchSpeakersDistMeta}.
}
\usage{data(dutchSpeakersDist)}
\format{
  A data frame for a 165 by 165 matrix of between-speaker differences.
}
\source{
  http://lands.let.kun.nl/cgn/
  data collected and analyzed in collaboration with Patrick Juola
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
    data(dutchSpeakersDist)
    dutchSpeakersDist.d = as.dist(dutchSpeakersDist)
    dutchSpeakersDist.mds = cmdscale(dutchSpeakersDist.d, k = 3)

    data(dutchSpeakersDistMeta)
    dat = data.frame(dutchSpeakersDist.mds, 
       Sex = dutchSpeakersDistMeta$Sex, 
       Year = dutchSpeakersDistMeta$AgeYear, 
       EduLevel = dutchSpeakersDistMeta$EduLevel)
    dat = dat[!is.na(dat$Year),]

    par(mfrow=c(1,2))
    plot(dat$Year, dat$X1, xlab="year of birth", 
       ylab = "dimension 1", type = "p")
    lines(lowess(dat$Year, dat$X1))
    boxplot(dat$X3 ~ dat$Sex, ylab = "dimension 3")
    par(mfrow=c(1,1))

    cor.test(dat$X1, dat$Year, method="sp")
    t.test(dat$X3~dat$Sex)
	}
}
\keyword{datasets}
