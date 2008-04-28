`plotAll.fnc` <-
function(reslist, sameYrange=TRUE, ylabel, intrName=NA, pos="end", ylimit=NA,
 addlines=FALSE, cexsize = 0.6, conditioningVals=NA, conditioningColors=1, conditioningLines=1) {

   if (length(conditioningColors)==1) conditioningColors = rep(1, 1000)   # 1000 should be more than enough
   if (length(conditioningLines)==1) conditioningLines = rep(1, 1000)   # 1000 should be more than enough

   if (sameYrange) {
     ylimit = getRange.fnc(reslist)
   }

   if (is.na(pos)) pstn = 2
   else {
     if (pos=="beg") pstn = 4
     else pstn = 2
   }  # we will use pstn (position) to ensure that the string is adjusted away from the margin

   for (i in 1:length(reslist)) {
     if ((!sameYrange) & (length(ylimit)==1)) {
       ylimit = getRange.fnc(reslist[[i]])
     }
     if (is.data.frame(reslist[[i]])) {
       lst = reslist[[i]]
       n = 1
     } else {
       lst = reslist[[i]][[1]]
       n = length(reslist[[i]])
     }
     if ("Levels" %in% colnames(lst)) {
       isfactor = TRUE
     } else {
       isfactor = FALSE
     }
     if (lst$Type[1]==FALSE) {
       plot(lst$X, lst$Y, ylim=ylimit, type="l", 
         col = conditioningColors[1],
         lty = conditioningLines[1],
         xlab=as.character(lst$Predictor[1]), ylab=ylabel)
       if ("lower" %in% colnames(lst)) {
         lines(lst$X, lst$lower, lty=2, col = conditioningColors[1])
         lines(lst$X, lst$upper, lty=2, col = conditioningColors[1])
       }
       if (n>1) {
         if (!is.na(pos)) {
           ps = getPos.fnc(lst$Y, pos)
           epsilon = (max(ylimit)-min(ylimit))/40
           text(lst$X[ps], lst$Y[ps]+epsilon, 
             as.character(lst$Interaction[1]), cex=cexsize, pos=pstn)  
         }
         mtext(intrName, side=4, line=1, cex=cexsize, adj=0)
       }
     } else {  # a factor
       d = max(lst$X)-min(lst$X)
       xlimit = c(min(lst$X)-0.1*d, max(lst$X)+0.1*d)
       if (addlines) {
         plot(lst$X, lst$Y, ylim=ylimit, type="b", pch=21, xlim=xlimit,
         xlab=as.character(lst$Predictor[1]), ylab=ylabel, xaxt="n")
       } else {
         plot(lst$X, lst$Y, ylim=ylimit, type="p", pch=21, xlim=xlimit,
         xlab=as.character(lst$Predictor[1]), ylab=ylabel, xaxt="n")
       }
       mtext(lst$Levels, at=lst$X, side=1, line=1, cex=cexsize)
       if (n > 1) {
         if (!is.na(pos) & !is.na(conditioningVals)) {
           ps = getPos.fnc(lst$Y, pos)
           epsilon = (max(ylimit)-min(ylimit))/40
           text(lst$X[ps], lst$Y[ps]+epsilon, 
             labels=as.character(conditioningVals[1]), cex=cexsize, pos=pstn)  
         }
       }

       if ("lower" %in% colnames(lst)) {
         points(lst$X, lst$lower, lty=2,  pch="-")
         points(lst$X, lst$upper, lty=2,  pch="-")
       }
     }
     if (n > 1) {
       for (j in 2:n) {
         lst = reslist[[i]][[j]]
         if (lst$Type[1]==FALSE) {
           lines(lst$X, lst$Y, ylim=ylimit, type="l",  
             col=conditioningColors[j], lty=conditioningLines[j])
           if ("lower" %in% colnames(lst)) {
             lines(lst$X, lst$lower, lty=2, col = conditioningColors[j])
             lines(lst$X, lst$upper, lty=2, col = conditioningColors[j])
           }
           if (!is.na(pos[1]) & !is.na(conditioningVals[1])) {
             ps = getPos.fnc(lst$Y, pos)
             epsilon = (max(ylimit)-min(ylimit))/40
             text(lst$X[ps], lst$Y[ps]+epsilon, 
               labels=as.character(lst$Interaction[1]), cex=cexsize, pos=pstn)  
           }
         } else {
           if (addlines) {
             lines(lst$X, lst$Y, ylim=ylimit, type="b", pch=21, 
               col=conditioningColors[j], lty=conditioningLines[j], 
               xlab=as.character(lst$Predictor[1]), ylab=ylabel)
           } else {
             points(lst$X, lst$Y, ylim=ylimit, type="p", pch=21, 
             xlab=as.character(lst$Predictor[1]), ylab=ylabel)
           }
           mtext(intrName, side=4, line=1, cex=cexsize, adj=0)
           if (!is.na(pos) & !is.na(conditioningVals)) {
             ps = getPos.fnc(lst$Y, pos)
             epsilon = (max(ylimit)-min(ylimit))/40
             text(lst$X[ps], lst$Y[ps]+epsilon, 
               labels=as.character(conditioningVals[j]), cex=cexsize, pos=pstn)  
           }
           if ("lower" %in% colnames(lst)) {
             points(lst$X, lst$lower, lty=2,  pch="-")
             points(lst$X, lst$upper, lty=2,  pch="-")
           }
         }
       }
     }
   }
}
