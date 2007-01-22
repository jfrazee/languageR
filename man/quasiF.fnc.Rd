\name{quasiF.fnc}
\alias{quasiF.fnc}
\title{Quasi-F test}
\description{
  The textbook Quasi-F test for a split-plot design with subjects, items,
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
\seealso{ ~~objects to See Also as \code{\link{quasiF.split.plot.fnc}}.}
\examples{
  data(splitplot)

  splitplot.lm = lm(RT ~ SOA + Item + Subject + 
     SOA:Subject + Item:Subject, data = splitplot)
  splitplot.aov = anova(splitplot.lm)

  quasiF.fnc(splitplot.aov["SOA","Mean Sq"], 
     splitplot.aov["Item:Subject", "Mean Sq"],
     splitplot.aov["SOA:Subject", "Mean Sq"], 
     splitplot.aov["Item", "Mean Sq"],
     splitplot.aov["SOA","Df"], 
     splitplot.aov["Item:Subject", "Df"],
     splitplot.aov["SOA:Subject", "Df"], 
     splitplot.aov["Item", "Df"])

  # much simpler is
  quasiF.split.plot.fnc(splitplot)$quasiF

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
