\name{collin.fnc}
\alias{collin.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate condition number with intercept included}
\description{
Calculates the condition number with the intercept included, following
Belsley, Kuh and Welsch (1980).
}
\usage{
collin.fnc(data, colvector)
}
\arguments{
  \item{data}{A data frame.}
  \item{colvector}{A vector with the column numbers in the data frame for which
  the collinearity is to be assessed. Only numeric predictors allowed.}
}
\value{
  A list with components
  \item{svd}{Singular value decomposition}
  \item{cindex}{Condition indices}
  \item{cnumber}{The condition number}
  \item{pi}{The phi matrix}
}
\references{
Belsley, D. A.  and  Kuh, E. and  Welsch, R. E. (1980) \emph{Regression
Diagnostics. Identifying Influential Data and Sources of Collinearity}, Wiley
Series in Probability and Mathematical Statistics, New York.
}
\author{F. J. Tweedie}
\seealso{\code{\link{kappa}}}
\examples{
data(english)
collin.fnc(english[english$AgeSubj=="young",], 7:29)$cnumber
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
