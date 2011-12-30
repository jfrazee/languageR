aovlmer.fnc <- function(object, mcmc, which, noMCMC = FALSE, ...) {
  # to avoid namespace conflicts, we redefine the function
  data(aovlmerfnc)
  return(aovlmer.fnc(object=object, mcmc=mcmc, which=which, noMCMC=noMCMC, ...))
}
