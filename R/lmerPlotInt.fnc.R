lmerPlotInt.fnc = function(lmermodel, xname, yname, intxyname, 
qntls=seq(0,1,by=0.1), view = 30, addStdError = FALSE, ndigits=2,
  nlev=30, which="matplot", shadow = 0.5, colour = "lightblue"){

  require("lme4", quietly = TRUE, character = TRUE)

  if (!inherits(lmermodel, "lmer")) stop("model object must be fitted with lmer")
  f <- function(x,y) {return(intercept + slopeX*x + slopeY*y + interactionxy*x*y)}

  coefs = fixef(lmermodel)
  intercept = coefs["(Intercept)"]

  if (xname %in% names(coefs)) slopeX = coefs[xname]
  else stop(paste(xname, "is not a valid predictor name"))

  if (yname %in% names(coefs)) slopeY = coefs[yname]
  else stop(paste(yname, "is not a valid predictor name"))

  if (intxyname %in% names(coefs)) interactionxy = coefs[intxyname]
  else stop(paste(intxyname, "is not a valid predictor name"))

  dat = lmermodel@frame

   
  x = unique(as.numeric(quantile(dat[,xname], qntls)))
  y = unique(as.numeric(quantile(dat[,yname], qntls)))
  z <- outer(x, y, f)

  if (addStdError) {
    sdError = lmermodel@devComp["scale"]
    z1 = matrix(rnorm(length(z), 0, sdError), nrow(z), ncol(z))
    z = z + z1
  }

  pars = par(c("mfrow", "mar"))
  if (which == "all") {
     par(mfrow=c(2,2), mar=c(4.5,4,1,1))
  }

  if ((which == "persp") | (which == "all")) {
    persp(x, y, z, theta=view, col=colour, ticktype = "detailed", 
           xlab=xname, ylab=yname,zlab=as.character(formula(lmermodel))[2],
           shade=shadow, phi=20)
  }     
  if ((which == "contour") | (which == "all")) {
    contour(x, y, z, col="blue", nlevel=nlev, xlab=xname, ylab=yname)
  } 
  if ((which == "matplot") | (which == "all")) {
    offset = (max(x)-min(x))/5
    matplot(x, z, xlab=xname,type="l",ylab=as.character(formula(lmermodel))[2],
       xlim=c(min(x), max(x)+offset))
    x1 = rep(x[length(x)],length(y))
    y1 = f(x1, y)
    if (which == "all") text(x1+offset, y1, round(y,ndigits), adj=+0.8, cex=0.8) 
    else  text(x1+offset, y1, round(y,ndigits), adj=+0.8, cex=0.7)
  }
  if ((which == "image") | (which == "all")) {
    image(x, y, z, col = heat.colors(10))
  }
  if (which == "all") par(pars)
}

