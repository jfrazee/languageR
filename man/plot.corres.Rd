\name{plot.corres}
\alias{plot.corres}
\title{Plot method for correspondence objects}
\description{
This function defines a plot method for correspondence objects.
}
\usage{
plot.corres(x, main = "", addcol = TRUE, extreme = 0, rcex = 1, rcol = 1, 
rlabels = "", stretch = 1.4, ccex = 1, ccol = 2, clabels = "", ...)
}
\arguments{
  \item{x}{A correspondence object as produced by \code{corres}.}
  \item{main}{A string used for the main title of the plot.}
  \item{addcol}{A logical, if true, columns are added to the plot.}
  \item{extreme}{If nonzero, defines quantiles that define the extremes
     such that only data points exceeding these extremes are plotted.}
  \item{rcex}{sets cex graphical parameter for rows.}
  \item{rcol}{sets color for rows.}
  \item{rlabels}{vector of row labels.}
  \item{stretch}{a number defining the degree to which the columns (or rows)
     should be stretched out for visual presentation.}
  \item{ccex}{sets cex graphical parameter for columns.}
  \item{ccol}{sets color for columns.}
  \item{clabels}{vector of column labels.}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  A plot on the graphics device.
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{corres.fnc}}, \code{link{corsup.fnc}}.}
\examples{
  data(oldFrench)
  data(oldFrenchMeta)

  oldFrench.ca = corres.fnc(oldFrench)

  plot(oldFrench.ca)

  plot(oldFrench.ca, rlabels = oldFrenchMeta$Genre, 
  rcol = as.numeric(oldFrenchMeta$Genre), rcex = 0.5, 
  extreme = 0.1, ccol = "blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
