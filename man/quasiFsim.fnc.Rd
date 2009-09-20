\name{quasiFsim.fnc}
\alias{quasiFsim.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quasi-F test for specific simple design}
\description{
  This function carries out a Quasi-F test for data with columns 
	labelled SOA, Subject, Item.  
  This function is called by simulate.quasif.fnc, and is
  not intended for general use.
}
\usage{
quasiFsim.fnc(dat)
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
    data(quasif)
    quasiFsim.fnc(quasif)$quasiF
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
