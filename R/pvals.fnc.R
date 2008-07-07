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

        mcmc = mcmcsamp(object, n = nsim)
        hpd = HPDinterval(mcmc)

        # here we make a table for the fixed effects

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

        
        # For the random effects table, we have to figure out the names. 
        # I assume that the order in which we have the random effects parameters
        # in mcmc$ST follows the order in the lmer model object@ST. 
        # We work through the successive list components, and for each list component,
        # we work through the matrix column by column, skipping the zero entries above
        # the main diagonal.

        isStdev = vector()
        # a boolean vector to distinguish between random intercepts and random slopes/contrasts
        # the names for the random effects are going to live in the vector nms
        theValues = vector()
        names(object@ST) = names(object@flist)
        for (i in 1:length(object@ST)) {
          m = object@ST[[i]]
          if (i == 1) {
            nms = rownames(m)
            nms = paste(rep(names(object@ST)[i], length(nms)), nms, sep=" ")
            for (j in 1:nrow(m)) {
              isStdev = c(isStdev, TRUE)
              theValues = c(theValues, m[j,1])
            }
          } else {
            tmp = rownames(m)
            tmp = paste(rep(names(object@ST)[i], length(tmp)), tmp, sep=" ")
            nms = c(nms, tmp)
            for (j in 1:nrow(m)) {
              isStdev = c(isStdev, TRUE)
              theValues = c(theValues, m[j,1])
            }
          }
          if (nrow(m) > 1) {
            for (j in 2:ncol(m)) {
              for (k in 1:nrow(m)) {
                if (m[k,j] != 0) {
                  nms = c(nms, paste(names(object@ST)[i], colnames(m)[j-1], rownames(m)[k], sep=" "))
                  isStdev = c(isStdev, FALSE)
                  theValues = c(theValues, m[k,j])
                }
              }
            }
          }
        }
        

        # we now combine the mean of the sampled values for sigma with the HPDs
        rowSigma = round(c(mean(mcmc@sigma), median(mcmc@sigma), as.vector(hpd$sigma)), ndigits)

        # backtransform to original scale (to be implemented)
        #for (i in 1:nrow(mcmc@ST)) 
        #  if (isStdev[i]) 
        #    mcmc@ST[i,] = mcmc@ST[i,]*mcmc@sigma

        mcmcHPD = as.data.frame(HPDinterval(mcmc)$ST)
        rownames(mcmcHPD) = nms
        means = apply(mcmc@ST,1,mean)
        medians = apply(mcmc@ST,1,median)

        # we can now combine means, medians and HPDs for the random effects on the original scales
        rowsST = round(cbind(means, medians, mcmcHPD), ndigits)
        rownames(rowsST) = nms
        random = data.frame(rbind(rowsST,rowSigma)) 
        rownames(random)[nrow(random)] = "Sigma"
        colnames(random) = c("MCMCmean", "MCMCmedian", "HPD95lower", "HPD95upper")
        #random$Estimate = c(theValues, sgma)
        #random = random[,c(5,1:4)]

        mcmcM = as.matrix(mcmc)

        if (addPlot) {
          beg = ncol(mcmcM)-length(nms)
          end = ncol(mcmcM)
          colnames(mcmcM)[beg:end] = c(nms, "sigma")
          m = data.frame(Value = mcmcM[,1], Predictor = rep(colnames(mcmcM)[1], nrow(mcmcM)))
          for (i in 2:ncol(mcmcM)) {
            mtmp = data.frame(Value = mcmcM[,i], Predictor = rep(colnames(mcmcM)[i], nrow(mcmcM)))
            #if (logScale) 
            #  if (i %in% (beg:end)) 
            #    if ((i == end) | (isStdev[i-beg+1])) 
            #      mtmp = data.frame(Value = 2*log(mcmcM[,i]), Predictor = rep(colnames(mcmcM)[i], nrow(mcmcM)))
            m = rbind(m, mtmp)
          }
          print(densityplot(~Value|Predictor, data=m, scales=list(relation="free"),
          par.strip.text = list(cex = 0.75), xlab="Posterior Values", ylab = "Density", pch="."))

          if (withMCMC) {
            return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
            random = format(random, digits=ndigits, sci=FALSE), mcmc=as.data.frame(mcmcM)))
          } else {
            return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
            random = format(random, digits=ndigits, sci=FALSE)))
          }
        } else {
          return(list(fixed = format(fixed, digits=ndigits, sci=FALSE), 
            random = format(random, digits=ndigits, sci=FALSE)))
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

      cat("the input model is not an lmer, glmer or mer object\n")
      return()

  }
}
