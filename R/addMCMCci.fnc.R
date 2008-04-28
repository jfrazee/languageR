`addMCMCci.fnc` <-
function(mcmcM, model, m, fun, pred, predname=NA, factor=FALSE) {
  if (is.matrix(mcmcM)) {
    # get mcmc-derived confidence intervals
    confints = getMCMCintervals.fnc(fixf=fixef(model), 
      mcmcMatrix=mcmcM, m=m)
    if (is.na(predname)) predname = pred
    if (!factor) {
      x = m[,predname]
    } else {
      x = 1:nrow(m)
    }
    confints$X=x
    dfr = data.frame(x = x,
          lower=transforming.fnc(confints[,1], fun), 
          upper = transforming.fnc(confints[,2], fun))
    return(dfr)
  } else {
    if (!is.na(mcmcM[1]))
      stop("warning: mcmcM argument to addMCMCci.fnc is not a matrix\n")
  }
}

