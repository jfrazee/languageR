`simulateQuasif.fnc` <-
function(dat,  with = TRUE, nruns = 100, nsub = NA, nitem = NA, ...) {
  require("MASS", quietly = TRUE, character = TRUE)
  require("lme4", keep.source = F, quietly = TRUE, character = TRUE)

  if (with) {
     model.lmer = lmer(RT ~ SOA + (1 + SOA|Subject) + (1|Item), 
                       data = dat)
  } else {
     model.lmer = lmer(RT ~ 1 + (1 + SOA|Subject) + (1|Item), 
                       data = dat)
  }
  res = matrix(0, nruns, 5)  # 5 p-values for each simulation run
  for (run in 1:nruns) {     # loop over simulation runs

    # create simulated reaction times using the fixed-effects estimates
    # and random effect parameters estimated from ranef() and resid()
    # for a model fitted to the input data set

    # changing the power of the experiment by increasing subjects and items 
    if (is.na(nsub)) nsub = nlevels(dat$Subject)
    if (is.na(nitem)) nitem = nlevels(dat$Item)

    # we extract the parameters and build vectors of random effects
    su = mvrnorm(nsub, mu = c(0, 0), 
         Sigma = cov(ranef(model.lmer)$Subject)) 
    subjects = data.frame(Subject = paste("S", 1:nsub, sep=""),
         rSubIntercept = su[,1], rSubSOA = su[, 2])
    items = data.frame(Item = paste("W", 1:nitem, sep=""),
         rItemIntercept = rnorm(nitem, 0, sd(ranef(model.lmer)$Item)))
    intercept = fixef(model.lmer)[1]
    if (with) slope = fixef(model.lmer)[2] else slope = 0
      
    # create a fully crossed design, after which we delete items from SOA conditions
    # to obtain the required nesting
    simdat = expand.grid(Subject = paste("S", 1:nsub, sep=""),
                Item = paste("W", 1:nitem, sep=""),
                SOA = c("short", "long"))
    simdat$SOA = relevel(simdat$SOA, ref="long")
    itemsvec = as.character(simdat$Item)
    out1 = which(as.numeric(substr(itemsvec,2,nchar(itemsvec))) <= nitem/2 & simdat$SOA == "long")
    out2 = which(as.numeric(substr(itemsvec,2,nchar(itemsvec))) > nitem/2 & simdat$SOA != "long")
    simdat = simdat[-c(out1, out2),]
    rownames(simdat) = 1:nrow(simdat)

    # and add fixed and random effects components
    simdat$population = intercept + as.numeric(simdat$SOA=="short")*slope
    simdat = merge(simdat, subjects[,1:3], by.x="Subject", by.y="Subject")
    simdat$rSubSOA =  as.numeric(simdat$SOA=="short") * simdat$rSubSOA
    simdat = merge(simdat, items[,1:2], by.x="Item", by.y="Item")
    simdat$error = rnorm(nrow(simdat), 0, sd(resid(model.lmer)))

    # calculate the simulated RTs
    simdat$RTsim = apply(simdat[,4:8], 1, sum)

    # finally, run the models
    sim.lmer = lmer(RTsim ~ SOA +  (1 + SOA|Subject) + (1|Item), 
                    data = simdat)
    pvalues = pvals.fnc(sim.lmer, nsim = 10000)$fixed
    res[run, 1] = as.numeric(pvalues["SOAshort",6])
    res[run, 2] = quasiFsim.fnc(simdat)$p
    res[run, 3] = subjects.quasif.fnc(simdat)$p
    res[run, 4] = items.quasif.fnc(simdat)$p
    res[run, 5] = as.numeric(pvalues["SOAshort", 5])
    cat(".")  # prints dot on command line - a simple progress report
  }
  cat("\n")   # reset cursor to the beginning of the next line

  res = data.frame(res)   # turn the matrix into a data frame
  colnames(res) = c("lmer", "quasi-F", "by-subject", "by-item", "pMCMC")
  
  # calculate proportion of cases where p-value is smaller than 0.05
  alpha05 = c(apply(res[,1:5] < 0.05, 2, sum)/nrow(res),
              sum(apply(res[,3:4] < 0.05, 1, sum) == 2)/nrow(res))
  names(alpha05)[6] = "F1+F2"

  # calculate proportion of cases where p-value is smaller than 0.01
  alpha01 = c(apply(res[, 1:5] < 0.01, 2, sum)/nrow(res),
              sum(apply(res[,3:4] < 0.01, 1, sum) == 2)/nrow(res))
  names(alpha01)[6] = "F1+F2"
  alpha01 = alpha01[c(2,3,4,6,1,5)]
  alpha05 = alpha05[c(2,3,4,6,1,5)]
	names(alpha01)[5] = "lmer:pt"
	names(alpha01)[6] = "lmer:pMCMC"
	names(alpha05)[5] = "lmer:pt"
	names(alpha05)[6] = "lmer:pMCMC"
  return(list(alpha05 = alpha05, alpha01 = alpha01, res = res, with = with, last = simdat))
}

