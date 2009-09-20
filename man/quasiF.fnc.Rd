\name{quasiF.fnc}
\alias{quasiF.fnc}
\title{Quasi-F test}
\description{
  The textbook Quasi-F test for a design with subjects, items,
  and a single factorial predictor.  Included for educational purposes
  for this specific design only. 
}
\usage{
quasiF.fnc(ms1, ms2, ms3, ms4, df1, df2, df3, df4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms1}{Mean squares Factor}
  \item{ms2}{Mean squares Item:Subject}
  \item{ms3}{Mean squares Factor:Subject}
  \item{ms4}{Mean squares Item}
  \item{df1}{Degrees of freedom Factor}
  \item{df2}{Degrees of freedom Item:Subject}
  \item{df3}{Degrees of freedom Factor:Subject}
  \item{df4}{Degrees of freedom Item}
}
\value{
  A list with components
  \item{F }{Quasi-F value.}
  \item{df1}{degrees of freedom numerator.}
  \item{df2}{degrees of freedom denominator.}
  \item{p}{p-value.}
}
\author{R. H. Baayen}
\seealso{ ~~objects to See Also as \code{\link{quasiFsim.fnc}}.}
\examples{
  data(quasif)

  quasif.lm = lm(RT ~ SOA + Item + Subject + 
     SOA:Subject + Item:Subject, data = quasif)
  quasif.aov = anova(quasif.lm)

  quasiF.fnc(quasif.aov["SOA","Mean Sq"], 
     quasif.aov["Item:Subject", "Mean Sq"],
     quasif.aov["SOA:Subject", "Mean Sq"], 
     quasif.aov["Item", "Mean Sq"],
     quasif.aov["SOA","Df"], 
     quasif.aov["Item:Subject", "Df"],
     quasif.aov["SOA:Subject", "Df"], 
     quasif.aov["Item", "Df"])

  # much simpler is
  quasiFsim.fnc(quasif)$quasiF

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
