pvals.fnc <- function(object, nsim=10000, ndigits = 4, withMCMC = FALSE, addPlot = TRUE, 
   ...) {

  require("lme4", quietly = TRUE, character = TRUE)

  # tested for  lme4 "lme4_0.999375-20.tar.gz" and for "R version 2.7.1 (2008-07-07)"


  if (is(object, "mer")) {  

      coefs = summary(object)@coefs
      ncoef = length(coefs[,1])
      sgma = summary(object)@sigma
      
      if (nsim > 0) {   # nsim is n for mcmcsamp

        if (colnames(coefs)[3] == "z value") {
          stop("mcmc sampling is not yet implemented for generalized mixed models\n")
          # the error message I get when I try to apply mcmcsamp() is:
          #   Error in .local(object, n, verbose, ...) : Update not yet written
        }

        mcmc = try(mcmcsamp(object, n=nsim),silent=TRUE)
        if (is(mcmc, "try-error")) {
          stop("MCMC sampling is not yet implemented in lme4_0.999375-20\n  for models with random correlation parameters\n")
          # stop("cannot carry out MCMC sampling:\nCode for non-trivial theta_T not yet written\n")
        }

        hpd = HPDinterval(mcmc)

        # the table for the fixed effects

        mcmcfixef = t(mcmc@fixef)
        nr <- nrow(mcmcfixef)
        prop <- colSums(mcmcfixef > 0)/nr
        ans <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))

        fixed = data.frame(
          Estimate   = round(as.numeric(coefs[,1]), ndigits), 
          MCMCmean   = round(apply(t(mcmc@fixef),2,mean), ndigits),
          HPD95lower = round(hpd$fixef[,1], ndigits),
          HPD95upper = round(hpd$fixef[,2], ndigits),
          pMCMC      = round(ans, ndigits), 
          pT         = round(2*(1-pt(abs(coefs[,3]),nrow(object@frame)-ncoef)), ndigits),
          row.names  = names(coefs[,1])
        )
        colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"

        # the table for the random effects

        ranefNames = names(object@flist)
        assigned = attr(object@flist, "assign")
        n = length(assigned)+1

        dfr = data.frame(
          Groups   = rep("", n),    
          Name     = rep("", n),
          Std.Dev. = rep(0, n),
          MCMCmedian   = rep(0, n),
          MCMCmean     = rep(0, n),
          HPD95lower    = rep(0, n),
          HPD95upper    = rep(0, n)
        )
        dfr$Groups = as.character(dfr$Groups)
        dfr$Name = as.character(dfr$Name)

        for (i in 1:length(object@ST)) {
          dfr$Groups[i] = ranefNames[assigned[i]]
          dfr$Name[i] = colnames(object@ST[[i]])
          dfr$Std.Dev.[i] = round(object@ST[[i]]*sgma, ndigits)
          dfr$MCMCmedian[i] = round(median(mcmc@ST[i,]*mcmc@sigma), ndigits)
          dfr$MCMCmean[i] = round(mean(mcmc@ST[i,]*mcmc@sigma), ndigits)
          hpdint = as.numeric(HPDinterval(mcmc@ST[i,]*mcmc@sigma))
          dfr$HPD95lower[i] = round(hpdint[1], ndigits)
          dfr$HPD95upper[i] = round(hpdint[2], ndigits)
        }
        dfr[n,1] = "Residual"
        dfr[n,2] = " "
        dfr[n,3] = round(sgma, ndigits)
        dfr[n,4] = round(median(mcmc@sigma), ndigits)
        dfr[n,5] = round(mean(mcmc@sigma), ndigits)
        hpdint = as.numeric(HPDinterval(mcmc@sigma))
        dfr[n,6] = round(hpdint[1], ndigits)
        dfr[n,7] = round(hpdint[2], ndigits)

        mcmcM = as.matrix(mcmc)
        k = 0
        for (j in (ncol(mcmcM)-n+1):(ncol(mcmcM)-1)) {
          k = k + 1
          mcmcM[,j] = mcmcM[,j]*mcmcM[,"sigma"]
          colnames(mcmcM)[j] = paste(dfr$Group[k], dfr$Name[k], sep=" ")
        }

        if (addPlot) {
          m = data.frame(Value = mcmcM[,1], Predictor = rep(colnames(mcmcM)[1], nrow(mcmcM)))
          for (i in 2:ncol(mcmcM)) {
            mtmp = data.frame(Value = mcmcM[,i], Predictor = rep(colnames(mcmcM)[i], nrow(mcmcM)))
            m = rbind(m, mtmp)
          }
          print(densityplot(~Value|Predictor, data=m, scales=list(relation="free"),
          par.strip.text = list(cex = 0.75), xlab="Posterior Values", ylab = "Density", pch="."))

          if (withMCMC) {
            return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
            random = dfr, mcmc=as.data.frame(mcmcM)))
          } else {
            return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
            random = dfr))
          }
        } else {
          return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
            random = dfr))
        }

      } else {

        coefs = summary(object)@coefs
        ncoef = length(coefs[,1])

        fixed = data.frame(
          Estimate   = round(as.numeric(coefs[,1]), ndigits), 
          pT         = round(2*(1-pt(abs(coefs[,3]),nrow(object@frame)-ncoef)), ndigits),
          row.names  = names(coefs[,1]))
          colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"

        return(list(fixed = format(fixed, digits=ndigits, sci=FALSE)))
      }

  } else {

      cat("the input model is not a mer object\n")
      return()

  }
}
