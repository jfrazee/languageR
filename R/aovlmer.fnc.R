aovlmer.fnc <- function(object, mcmc, which, noMCMC = FALSE, ...) {

  require("lme4", quietly = TRUE, character = TRUE)

  if (!(is(object, "mer"))) 
		stop("first argument should be a mer model object")


  sumry = summary(object)
  nobs  = nrow(object@frame)
  ncoef = nrow(sumry@coefs)
 
  #options(show.error.messages=F)
  #anov = try(anova(object))
  #options(show.error.messages=T)
  #if (class(anov) != "try-error") {
  #  anov$F = anov[,"Mean Sq"]/attr(sumry,"sigma")^2
  #  anov$Df2 = nobs-ncoef
  #  anov$p = 1-pf(anov$F, anov$Df, anov$Df2)
  #} else {
  #  anov = "anova function failed"
  #}

	if (noMCMC) {
    stop("single-argument anova() no longer supported by lmer()\nuse sequential likelihood ratio tests, or use MCMC sampling\n")
	} else {
    if (!is(mcmc, "data.frame")) stop("second argument should be a data frame")
    mcmc = mcmc[, which]      # next lines by Douglas Bates
    std <- backsolve(chol(var(mcmc)),
                     cbind(0, t(mcmc)) - colMeans(mcmc),
                     transpose = TRUE)
    sqdist <- colSums(std * std)
    pmcmc = sum(sqdist[-1] > sqdist[1])/nrow(mcmc)

    return(list(p=pmcmc, which = which))
  }
}
