aovlmer.fnc <- function(object, mcmc, which, noMCMC = FALSE, ...) {

  if (!is(object, "lmer")) stop("first argument should be an lmer model object")

  sumry = summary(object)
  nobs  = nrow(object@frame)
  ncoef = nrow(sumry@coefs)
 
  anov   = anova(object)
  anov$F = anov[,"Mean Sq"]/attr(sumry,"sigma")^2
  anov$Df2 = nobs-ncoef
  anov$p = 1-pf(anov$F, anov$Df, anov$Df2)

	if (noMCMC) {
		return(anov)
	} else {
    if (!is(mcmc, "mcmc")) stop("second argument should be an mcmc object")

    mcmc = mcmc[, which]      # next lines by Douglas Bates
    std <- backsolve(chol(var(mcmc)),
                     cbind(0, t(mcmc)) - colMeans(mcmc),
                     transpose = TRUE)
    sqdist <- colSums(std * std)
    pmcmc = sum(sqdist[-1] > sqdist[1])/nrow(mcmc)

    return(list(MCMC = list(p=pmcmc, which = which), Ftests = anov))
  }
}
