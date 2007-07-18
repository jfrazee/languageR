\name{lmerPlotInt.fnc}
\alias{lmerPlotInt.fnc}
\title{Plot the interaction of two linear numeric predictors in a model 
fitted with lmer}
\description{
 Visualization of an interaction in a model fitted with lmer of two 
 numeric predictors.
}
\usage{
lmerPlotInt.fnc(lmermodel, xname, yname, intxyname, qntls = seq(0, 1, by = 0.1), view = 30, addStdError = FALSE, ndigits = 2, nlev = 30, which = "matplot", shadow = 0.5, colour = "lightblue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmermodel}{an lmer model object}
  \item{xname}{name (character string) of first numeric predictor}
  \item{yname}{name (character string) of second numeric predictor}
  \item{intxyname}{name (character string) of the interaction in the lmer summary}
  \item{qntls}{vector of values to be shown for the second numeric predictor, defaults to deciles}
  \item{view}{specifies the viewing parameter \code{theta} for the perspective plot}
  \item{addStdError}{add noise with the standard deviation of the residual error in the lmer model to the plot}
  \item{ndigits}{number of digits to show for the second numeric predictor }
  \item{nlev}{number of levels for the contour plot}
  \item{which}{choices are "matplot" (default), "contour", "persp", "image", 
	 and "all", in which case a 2 by 2 panel is shown with all four plots}
  \item{shadow}{the amount of \code{shade} for the perspective plot}
  \item{colour}{the color used for the perspective plot, defaults to "lightblue"}
}
\value{
  A plot is shown on the graphics device.
}
\author{R. H. Baayen}
\section{Warning}{
	This function should not be used to plot interactions when one of the predictors also has quadratic or higher terms in the model.
}
\examples{
	lexdec.lmer = lmer(RT~BNCw*Frequency+(1|Subject)+(1|Word), data=lexdec)
	lmerPlotInt.fnc(lexdec.lmer, "BNCw", "Frequency", "BNCw:Frequency", which="matplot")
}
\keyword{hplot}
