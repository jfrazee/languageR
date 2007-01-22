`plot.logistic.fit.fnc` <-
function(x = dative.glmm, data = dative, ...) {
  require(lme4, quietly = TRUE)
  require(Design, quietly = TRUE)
  if (class(x)[1]=="glmer") {
    depvar = as.character(formula(attr(x@frame, "terms")))[2]
    probs = 1/(1+exp(-fitted(x)))
  } else {
    if (class(x)[1] == "lrm") {
      depvar = as.character(formula(x$call))[2]
      probs = predict(x,type="fitted")
    } else {
      stop("first argument is not an lmer or lrm model")
    }
  }
  classes = cut2(probs, seq(0, 1, by = 0.1), levels.mean = TRUE)
  means = tapply(as.numeric(data[,depvar])-1, classes, mean)
  plot(as.numeric(names(means)), means, 
     xlab = "mean predicted probabilities", 
     ylab = "observed proportions")
  abline(0, 1, col = "grey")
  mtext(paste("R-squared: ", round(cor(as.numeric(names(means)), means)^2,2),
  sep=""), 3, 1.5)
}

