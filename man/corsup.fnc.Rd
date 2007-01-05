\name{corsup.fnc}
\alias{corsup.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Supplementary rows or columns in correspondence analysis}
\description{
Corsup calculates supplementary rows or columns for correspondence analysis.
}
\usage{
corsup.fnc(corres, bycol = TRUE, supp, plot = TRUE, font = 3, labels = "", 
cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corres}{A correspondence object.}
  \item{bycol}{A logical value indicating whether supplementary columns
       (the default) or supplementary rows are required.}
  \item{supp}{Supplementary rows or columns from a data frame with the
        same structure as the data frame used for the \code{corres} object.}
  \item{plot}{A logical value indicating whether supplementary rows or columns
        should be added to an already existing plot.}
  \item{font}{An integer specifying the font to be used for plotting.}
  \item{labels}{A character vector with row or column names to be used for plotting.}
  \item{cex}{A real specifying the font size required for plotting.}
}
\value{
  If \code{plot = FALSE}, a matrix with the supplementary coordinates.
  Otherwise, supplementary rows or columns are added to an already
  existing plot of a correspondence object.
}
\references{
F. Murtagh (2005) \emph{Correspondence Analysis and Data Coding with 
JAVA and R}, Chapman & Hall/CRC, Boca Raton, FL.
}
\author{
Extension of the code in Murtagh (2005) by R. Harald Baayen

Radboud University Nijmegen & Max Planck Institute for Psycholinguistics

Nijmegen, The Netherlands

email: baayen@mpi.nl
}
\seealso{\code{\link{corres.fnc}}}
\examples{
 \dontrun{
  data(oldFrench)
  data(oldFrenchMeta)
  prose = oldFrench[oldFrenchMeta$Genre=="prose" & 
          !is.na(oldFrenchMeta$Year),]
  proseinfo = oldFrenchMeta[oldFrenchMeta$Genre=="prose" & 
          !is.na(oldFrenchMeta$Year),]
  proseinfo$Period = as.factor(proseinfo$Year <= 1250)

  prose.ca = corres.fnc(prose)
  plot(prose.ca, addcol = FALSE, 
  rcol = as.numeric(proseinfo$Period) + 1, 
  rlabels = proseinfo$Year, rcex = 0.7)

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
