\name{plotAll.fnc}
\alias{plotAll.fnc}
\title{create plot or plots for list with data frames for plot or subplots }
\description{
  given a list with one or more data frames with values for a plot (or subplot),
  create the actual plots
}
\usage{
plotAll.fnc(reslist, sameYrange = TRUE, ylabel, xlabel = NA, intrName = NA, 
  pos = "end", ylimit = NA, addlines=FALSE, cexsize = 0.6, conditioningVals=NA,
  conditioningColors=1, conditioningLines=1, lineColor=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reslist}{list with as components either a data frame or a list with data frames, the data frames specify X and Y coordinates and HPD intervals}
  \item{sameYrange}{ logical, if \code{TRUE}, the \code{ylim} for each panel will be chosen to accomodate the range of values across all panels in the plot }
  \item{ylabel}{ label to be used for the vertical axis}
  \item{xlabel}{ label to be used for the horizontal axis; this option is active
    only when a single predictor is plotted}
  \item{intrName}{ label for the interaction predictor, if present }
  \item{pos}{ location of legend values for interaction }
  \item{ylimit}{ if specified, overrides \code{sameYrange} for \code{ylim}}
  \item{addlines}{ if TRUE, adds line between levels of same factor(s)}
  \item{cexsize}{character expansion size for information in the plot for
    interactions, default is 0.6}
  \item{conditioningVals}{ vector of names of the levels of the conditioning
    factor in the interaction (the factor with different lines in the plot) }
  \item{conditioningColors}{ vector of names of the colors to be used
    for the levels of the conditioning factor in the interaction (the 
    factor with different lines in the plot) }
  \item{conditioningLines}{ vector of names of the line types to be used
    for the levels of the conditioning factor in the interaction (the 
    factor with different lines in the plot), by default solid lines }
  \item{lineColor}{ name of color to be used for the lines in the plot}
  \item{...}{ further graphical parameters to be passed down, none are currently implemented}
}
\details{
  Note that \code{reslist} may contain as elements lists of data frames, these then specify the separate points or lines to be plotted for a given interaction
}
\value{
  A plot is produced on the graphics device.
}
\references{ ~put references to the literature/web site here ~ }
\author{R. H. Baayen}
\note{
 not intended for independent use
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
  \dontrun{not intended for independent use}
}
\keyword{ regression }
