`plot.logistic.fit.fnc` <-
function(x, data, method = "cut", where = seq(0, 1, by=0.1), ...) {
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
  if (method == "cut") {
    classes = cut2(probs, where, levels.mean = TRUE)
    means = tapply(as.numeric(data[,depvar])-1, classes, mean)
  } else {
    if (method == "shingle") {
      sh = equal.count(probs)
      means = rep(0, length(levels(sh)))
      midpoints = rep(0, length(means))
      for (i in 1:length(levels(sh))) {
        means[i] = mean(probs[probs>levels(sh)[[i]][1] & probs < levels(sh)[[i]][2]])
        midpoints[i] = as.character(mean(levels(sh)[[i]]))
      }
      names(means) = as.character(midpoints)
    }
  }
  plot(as.numeric(names(means)), means, 
     xlab = "mean predicted probabilities", 
     ylab = "observed proportions", type="n")
  abline(0, 1, col = "grey")
  points(as.numeric(names(means)), means)
  mtext(paste("R-squared: ", round(cor(as.numeric(names(means)), means)^2,2),
  sep=""), 3, 1.5)
}

