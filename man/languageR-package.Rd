\name{languageR-package}
\alias{languageR-package}
\alias{languageR}
\docType{package}
\title{
Data sets and functions for 'Analyzing Linguistic Data'
}
\description{
Data sets and functions accompanying 'Analyzing Linguistic Data:
A practical introduction to statistics', Cambridge University Press, 2007.
}
\details{
\tabular{ll}{
Package: \tab languageR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-01-15\cr
License: \tab GNU public license \cr
}
The main function of this package is to make available the data sets
discussed and analyzed in 'Analyzing Linguistic Data:
A practical introduction to statistics using R', to appear with 
Cambridge University Press.  The following packages should be
installed, as ancillary functions in this package depend on them.
\describe{
    \item{\code{zipfR}}{for word frequency distributions}
    \item{\code{lme4}}{for mixed-effects models}
    \item{\code{coda}}{for Markov-Chain Monte Carlo estimation}
    \item{\code{lattice}}{for trellis graphics}
    \item{\code{Matrix}}{for mixed-effects modeling}
}
The following packages need to be installed for working through specific 
examples.
\describe{
    \item{\code{rms}}{for regression modeling}
    \item{\code{rpart}}{for CART trees}
    \item{\code{e1071}}{for support vector machines}
    \item{\code{MASS}}{for many useful functions}
    \item{\code{ape}}{for phylogenetic clustering}
}

The main convenience functions in this library are, by category:

\describe{
    \item{correspondence analysis}{(extending code by Murtagh, 2005)}
       \describe{
          \item{\code{corres.fnc}}{correspondence analysis}
          \item{\code{corsup.fnc}}{supplementary data}
       }
    \item{vocabulary richness}{(supplementing current zipfR functionality)}
       \describe{
          \item{\code{compare.richness.fnc}}{for two texts, compare richness}
          \item{\code{growth.fnc}}{empirical vocabulary growth data for text}
          \item{\code{growth2vgc}}{conversion to vgc object of zipfR}
          \item{\code{spectrum.fnc}}{creates frequency spectrum}
          \item{\code{text2spc.fnc}}{conversion to spc object of zipfR}
       }
    \item{lmer functions}{(p-values for mixed-effects models with lme4)}
		   \describe{
					\item{\code{pvals.fnc}}{p-values for table of coefficients including
					      MCMC}
					\item{\code{aovlmer.fnc}}{p-values for anova tables and/or MCMC
					      p-value for specified factor}
			 }
    \item{simulation functions}{(for comparing mixed models with traditional techniques including F1, F2, and F1+F2)}
       \describe{
           \item{\code{simulateRegression.fnc}}{simulate simple regression design}
           \item{\code{simulateQuasif.fnc}}{simulate data for Quasi-F ratios}
           \item{\code{simulateLatinsquare.fnc}}{simulating simple Latin-square design} 
       }
    \item{miscellaneous}{(convenience functions)}
       \describe{
          \item{\code{pairscor.fnc}}{scatterplot matrix with correlation tests}
          \item{\code{collin.fnc}}{collinearity diagnostics}
          \item{\code{pvals.fnc}}{p-values and MCMC confidence intervals for mixed models}
          \item{\code{plot.logistic.fit.fnc}}{diagnostic visualization for logistic models}
          \item{\code{xylowess.fnc}}{trellis scatterplots with smoother}
          \item{\code{mvrnormplot.fnc}}{scatterplot for bivariate standard normal random numbers with regression line}
					\item{\code{lmerPlotInt.fnc}}{offers choice of four ways to visualize an interaction between two numeric predictors in an lmer model}
       }
}
}
\author{
R. H. Baayen

University of Alberta, Edmonton, Canada

\email{harald.baayen@gmail.com}

Maintainer:  harald.baayen@gmail.com 
}
\references{
R. H. Baayen (2007) \emph{Analyzing Linguistic Data: A practical introduction
to statistics using R}, Cambridge: Cambridge University Press. 
}
\keyword{package}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{\dontrun{
  library(languageR)
  data(package="languageR")
}}
