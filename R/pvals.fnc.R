pvals.fnc <- function(object, nsim=10000, ndigits = 4, withMCMC = FALSE, ...) {

  require("coda", quietly = TRUE, character = TRUE)
  require("lme4", quietly = TRUE, character = TRUE)

  if (!is(object, "lmer")) stop("argument should be an lmer model object")

  coefs = summary(object)@coefs
  ncoef = length(coefs[,1])

  if ((nsim > 0) & (class(object) == "lmer")) {
    mcmc = mcmcsamp(object, n = nsim)
    hpd = HPDinterval(mcmc)
    sumry = summary(mcmc)$statistics

    nr <- nrow(mcmc)
    prop <- colSums(mcmc[ , 1:ncoef] > 0)/nr
    ans <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))

    fixed = data.frame(
      Estimate   = round(as.numeric(coefs[,1]), ndigits), 
      MCMCmean   = round(sumry[,1][1:ncoef], ndigits),
      HPD95lower = round(hpd[1:ncoef,1], ndigits),
      HPD95upper = round(hpd[1:ncoef,2], ndigits),
      pMCMC      = round(ans, ndigits), 
      pT         = round(2*(1-pt(abs(coefs[,3]),nrow(object@frame)-2)), ndigits),
      row.names=names(coefs[,1])
    )
    colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
  
    v = (ncoef+1):nrow(hpd)
    random = data.frame(
      MCMCmean   = sumry[,1][v],
      HPD95lower = hpd[v,1],
      HPD95upper = hpd[v,2]
    )
  
    nms = rownames(random)
  
    logs = substr(rownames(random),1,3)
    rows = logs == "log"
    random[rows, ] = sqrt(exp(random[rows,]))
    nms[rows] = substr(nms[rows],5,nchar(nms)-1)
    nms[1] = substr(nms[1],1,nchar(nms[1])-2)  # strip ^2 from sigma^2
  
    atanhs = substr(rownames(random), 1, 5)
    rows = atanhs == "atanh"          
    if (sum(rows) > 0) {
      random[rows, ] = tanh(random[rows, ])
      nms[rows] = substr(nms[rows],7,nchar(nms[rows])-1)
    }
    rownames(random) = nms
  
    if (withMCMC)
      return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
      random = format(random, digits=ndigits, sci=FALSE), mcmc=mcmc))
    else
      return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
      random = format(random, digits=ndigits, sci=FALSE)))
  } else {
    fixed = data.frame(
      Estimate   = as.numeric(coefs[,1]), 
      pT         = round(2*(1-pt(abs(coefs[,3]),nrow(object@frame)-2)), ndigits),
      row.names=names(coefs[,1])
    )
    colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
    return(fixed)
  }
}
