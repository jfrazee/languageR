`getMCMCintervals.fnc` <-
function(fixf, mcmcMatrix, m) {
  nfixed = length(fixf)
  nsamples = nrow(mcmcMatrix)
  mat = matrix(0, nrow(m), nsamples)
  for (i in 1:nsamples) {
    mat[,i] = m %*% as.numeric(mcmcMatrix[i,1:nfixed])
  }
  # require("coda", quietly=TRUE)
  hpd = as.data.frame(HPDinterval(as.mcmc(t(mat))))
  # for lme4 1.0 probably:
  # hpd = as.data.frame(HPDinterval(t(mat)))
  return(hpd)
}

