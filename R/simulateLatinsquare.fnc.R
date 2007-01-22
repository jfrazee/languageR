`simulateLatinsquare.fnc` <-
function(dat, with = TRUE, trial = 0, nruns = 100, nsub = NA, nitem = NA, ...) {
  # dat: data frame that will be modelled with lmer() to obtain fixed and random
  #      effects, calculated from fixed() and ranef()
  # with: if TRUE then SOA effect is implemented, otherwise all contrasts zero
  # trial: if nonzero, for each successive trial, 'trial' ms will be added to the
  #      reaction time.  A sensible value is something lik 0.2.  Each subject
  #      responds to the trials in a different random order
  # nruns: number of simulation runs
  # nsub: if NA, number of subjects is determined using the subjects in dat
  #       if set by hand, should be a multiple of 3
  # nitem: if NA, number of items is determined using the items in dat
  #       if set by hand, should be a multiple of 3
  
  require("lme4", quietly = TRUE)
  require("coda", quietly = TRUE)


  if (with) {
     model.lmer = lmer(RT ~ SOA + (1|Subject) + (1|Word), 
                       data = dat)
  } else {
     model.lmer = lmer(RT ~ 1 + (1|Subject) + (1|Word), 
                       data = dat)
  }
  res = matrix(0, nruns, 3)  # 3 p-values for each simulation run

  # create simulated reaction times using the fixed-effects estimates
  # and random effect parameters estimated from ranef() and resid()
  # for a model fitted to the input data set

  # changing the power of the experiment by increasing subjects and items 
  if (is.na(nsub)) {
    nsub = nlevels(dat$Subject)
  } else {
    if (nsub %% 3 != 0) {
      cat("number of subjects should be multiple of 3\n")
      return(-1)
    }
  }
  if (is.na(nitem)) {
    nitem = nlevels(dat$Word)
    if (nitem %% 3 != 0) {
      cat("number of items should be multiple of 3\n")
      return(-1)
    }
  }

  
  for (run in 1:nruns) {     # loop over simulation runs

    # we extract the parameters and build vectors of random effects
    subjects = data.frame(Subject = paste("S", 1:nsub, sep=""),
       rSubIntercept = rnorm(nsub, 0, sd(ranef(model.lmer)$Subject)))
    items = data.frame(Word = paste("W", 1:nitem, sep=""),
       rItemIntercept = rnorm(nitem, 0, sd(ranef(model.lmer)$Word)))

    # create a fully crossed design, after which we delete items from SOA conditions
    # to obtain the required nesting
    simdat = expand.grid(Subject = paste("S", 1:nsub, sep=""),
                         Word    = paste("W", 1:nitem, sep=""),
                         SOA     = c("short", "medium", "long"))
    simdat$SOA = relevel(simdat$SOA, ref="long")

    # and now create Latin Square 
    itemsvec = as.character(simdat$Word)
    simdat$ItemCnt = as.numeric(substr(itemsvec, 2, nchar(itemsvec)))
    subjectsvec = as.character(simdat$Subject)
    simdat$SubjectCnt = as.numeric(substr(subjectsvec, 2, nchar(subjectsvec)))

    g1 = simdat[simdat$SOA == "short" & 
                simdat$ItemCnt %in% 1:(nitem/3) & 
                simdat$SubjectCnt %in% 1:(nsub/3),]
    g1$List = rep("L1", nrow(g1))
    g2 = simdat[simdat$SOA == "medium" & 
                simdat$ItemCnt %in% ((nitem/3)+1):(2*nitem/3) & 
                simdat$SubjectCnt %in% 1:(nsub/3),]
    g2$List = rep("L2", nrow(g1))
    g3 = simdat[simdat$SOA == "long" & 
                simdat$ItemCnt %in% ((2*nitem/3)+1):(3*nitem/3) & 
                simdat$SubjectCnt %in% 1:(nsub/3),]
    g3$List = rep("L3", nrow(g3))

    g4 = simdat[simdat$SOA == "medium" & 
                simdat$ItemCnt %in% 1:(nitem/3) & 
                simdat$SubjectCnt %in% ((nsub/3)+1):(2*nsub/3),]
    g4$List = rep("L1", nrow(g4))
    g5 = simdat[simdat$SOA == "long" & 
                simdat$ItemCnt %in% ((nitem/3)+1):(2*nitem/3) & 
                simdat$SubjectCnt %in% ((nsub/3)+1):(2*nsub/3),]
    g5$List = rep("L2", nrow(g5))
    g6 = simdat[simdat$SOA == "short" & 
                simdat$ItemCnt %in% ((2*nitem/3)+1):(3*nitem/3) & 
                simdat$SubjectCnt %in% ((nsub/3)+1):(2*nsub/3),]
    g6$List = rep("L3", nrow(g6))

    g7 = simdat[simdat$SOA == "long" & 
                simdat$ItemCnt %in% 1:(nitem/3) & 
                simdat$SubjectCnt %in% ((2*nsub/3)+1):(3*nsub/3),]
    g7$List = rep("L1", nrow(g1))
    g8 = simdat[simdat$SOA == "short" & 
                simdat$ItemCnt %in% ((nitem/3)+1):(2*nitem/3) & 
                simdat$SubjectCnt %in% (2*(nsub/3)+1):(3*nsub/3),]
    g8$List = rep("L2", nrow(g8))
    g9 = simdat[simdat$SOA == "medium" & 
                simdat$ItemCnt %in% ((2*nitem/3)+1):(3*nitem/3) & 
                simdat$SubjectCnt %in% ((2*nsub/3)+1):(3*nsub/3),]
    g9$List = rep("L3", nrow(g9))

    simdat=rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9)
    simdat$List = as.factor(simdat$List)
    rownames(simdat) = 1:nrow(simdat)
    simdat$Group = 
      as.factor(paste(rep("G", nrow(simdat)), rep(1:3, rep(nrow(simdat)/3,3)), sep=""))
    simdat = simdat[order(simdat$Subject),]
    n = nrow(simdat)/nsub
    trials = sample(1:n)
    for (i in 2:nsub) {
        trials = c(trials, sample(1:n))
    }
    simdat$Trial = trials

    # add simulated RTs to data frame

    # fixed effects 
    simdat$population = fixef(model.lmer)["(Intercept)"]

    if (with) {
       simdat[simdat$SOA=="medium",]$population = 
         simdat[simdat$SOA=="medium",]$population + fixef(model.lmer)["SOAmedium"]
       simdat[simdat$SOA=="short",]$population = 
         simdat[simdat$SOA=="short",]$population + fixef(model.lmer)["SOAshort"]
    }

    # and add random effects 
    simdat = merge(simdat, subjects, by.x="Subject", by.y="Subject")
    simdat = merge(simdat, items, by.x="Word", by.y="Word")
    simdat$error = rnorm(nrow(simdat), 0, sd(resid(model.lmer)))

    simdat$RTsim = apply(simdat[,9:12], 1, sum)
    if (trial > 0) {
      simdat$RTsim = simdat$RTsim + trial*simdat$Trial
      sim.lmer = lmer(RTsim ~ Trial + SOA + (1|Subject) + (1|Word), 
                 data = simdat)
    } else {
      sim.lmer = lmer(RTsim ~ SOA +  (1|Subject) + (1|Word), 
                 data = simdat)
    } 

    mcmc = mcmcsamp(sim.lmer, n = 10000)
		aov = aovlmer.fnc(sim.lmer, mcmc, c("SOAmedium", "SOAshort"))
		if (trial == 0) res[run, 1] = aov$Ftests$p
		else res[run,1] = aov$Ftests$p[2]
    res[run, 2] = aov$MCMC$p
    res[run, 3] = subjects.latinsquare.fnc(simdat)$p
    cat(".")  # prints dot on command line - a simple progress report
  }
  cat("\n")   # reset cursor to the beginning of the next line
  res = data.frame(res)   # turn the matrix into a data frame
  colnames(res) = c("Ftest", "MCMC", "F1")
  
  alpha05 = apply(res< 0.05, 2, sum)/nrow(res)
  alpha01 = apply(res< 0.01, 2, sum)/nrow(res)

  return(list(alpha05 = alpha05, alpha01 = alpha01, res = res, with = with))
}

