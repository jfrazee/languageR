\name{quasiF.split.plot.fnc}
\alias{quasiF.split.plot.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quasi-F test for simple split-plot design}
\description{
  This function carries out a Quasi-F test for a simple split-plot
  design for data with columns labelled SOA, Subject, Item.  
  This function is called by simulate.split.plot.fnc, and is
  not intended for general use.
}
\usage{
quasiF.split.plot.fnc(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame with RT (or RTsim), SOA, Subject and Item as
  predictors.}
}
\value{
  A list with components
  \item{p}{The p-value of the quasi-F test.}
  \item{data}{The input data.}
  \item{model}{The linear model fitted to the data.}
  \item{qF}{a list with F, df1, df2 and p-value of quasi-F test.}
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{quasiF.fnc}}.}
\examples{
    data(splitplot)
    quasiF.split.plot.fnc(splitplot)$quasiF
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
