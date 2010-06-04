\name{corres.fnc}
\alias{corres.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Correspondence Analysis}
\description{
Correspondence analysis for a contingency table. 
}
\usage{
corres.fnc(xtab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xtab}{A data frame cross-tabulating frequencies.}
}
\value{
  A correspondence object with summary and plot methods.
  The summary method lists eigenvalue rates and coordinates, correlations
  and contributions for Factor 1 and Factor 2.  By default, only the
  first six rows of the factor tables are shown.  Full tables are
  obtained by specifying \code{header = FALSE} when calling \code{summary}.
  For information on higher dimensions, set the option \code{n} to the
  desired number (e.g., \code{n = 3}) within \code{summary}.
  See \code{plot.corres} for documentation of plot options.
}
\references{
F. Murtagh (2005) \emph{Correspondence Analysis and Data Coding with 
JAVA and R}, Chapman & Hall/CRC, Boca Raton, FL.
}
\author{Extension of the code in Murtagh (2005) by R. Harald Baayen

Radboud University Nijmegen & Max Planck Institute for Psycholinguistics

Nijmegen, The Netherlands

email: baayen@mpi.nl
}
\seealso{See \code{\link{corsup.fnc}} for adding supplementary data to a 
correspondence plot, and \code{\link{plot.corres}} for plot options.}
\examples{
 \dontrun{
  data(oldFrench)
  oldFrench.ca = corres.fnc(oldFrench)
  oldFrench.ca
  summary(oldFrench.ca, head = TRUE)
  plot(oldFrench.ca)

  # more readable plot
  data(oldFrenchMeta)
  plot(oldFrench.ca, rlabels = oldFrenchMeta$Genre, 
  rcol = as.numeric(oldFrenchMeta$Genre), rcex = 0.5, 
  extreme = 0.1, ccol = "blue")

  # create subset of proze texts

  prose = oldFrench[oldFrenchMeta$Genre=="prose" & 
          !is.na(oldFrenchMeta$Year),]
  proseinfo = oldFrenchMeta[oldFrenchMeta$Genre=="prose" & 
          !is.na(oldFrenchMeta$Year),]
  proseinfo$Period = as.factor(proseinfo$Year <= 1250)

  prose.ca = corres.fnc(prose)
  plot(prose.ca, addcol = FALSE, 
  rcol = as.numeric(proseinfo$Period) + 1, 
  rlabels = proseinfo$Year, rcex = 0.7)

  # and add supplementary data for texts with unknown date of composition
  proseSup = oldFrench[oldFrenchMeta$Genre == "prose" & 
          is.na(oldFrenchMeta$Year),]
  corsup.fnc(prose.ca, bycol = FALSE, supp = proseSup, font = 2, 
  cex = 0.8, labels = substr(rownames(proseSup), 1, 4)) 
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
