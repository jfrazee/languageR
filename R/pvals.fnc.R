pvals.fnc <- function(object, nsim=10000, ndigits = 4, withMCMC = FALSE, addPlot = TRUE, ...) {
  # to avoid namespace conflicts, we redefine the function
  data(pvalsfnc)
  return(pvals.fnc(object=object, nsim=nsim, ndigits=ndigits, withMCMC=withMCMC, addPlot=addPlot))
}
