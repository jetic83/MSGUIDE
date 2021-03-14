
if (length(commandArgs()) > 4) {
  presSet <- as.numeric(commandArgs()[5])
} else {
  if (!exists("presSet")) presSet <- 1
}


predNamesSets = c(list(c("FINC", "VTNC", "BX_GS")))

avg.surv = function(survs) {
  avgsurv = survs[[1]]
  
  times_1 = c()
  times_2 = c()
  n = c(0, 0)
  for (i in 1:length(survs)) {
    n = n + survs[[i]]$n
    
    inds = c(1, which((survs[[i]]$time - survs[[i]]$time[-1])>0)[1])
    
    inds_1 = 1:(inds[2])
    inds_2 = min(inds[2]+1, length(survs[[i]]$time)):length(survs[[i]]$time)
    
    times_1 = c(times_1, survs[[i]]$time[inds_1])
    times_2 = c(times_2, survs[[i]]$time[inds_2])
  }
  
  times_1 = sort(unique(times_1))
  times_2 = sort(unique(times_2))
  
  events_1 = rep(0, length(times_1))
  events_2 = rep(0, length(times_2))
  
  censors_1 = rep(0, length(times_1))
  censors_2 = rep(0, length(times_2))
  
  for (i in 1:length(survs)) {
    inds = c(1, which((survs[[i]]$time - survs[[i]]$time[-1])>0)[1])
    
    inds_1 = 1:(inds[2])
    inds_2 = min(inds[2]+1, length(survs[[i]]$time)):length(survs[[i]]$time)
    
    ind = match(survs[[i]]$time[inds_1], times_1)
    events_1[ind] = events_1[ind] + survs[[i]]$n.event[inds_1]
    censors_1[ind] = censors_1[ind] + survs[[i]]$n.censor[inds_1]
    
    ind = match(survs[[i]]$time[inds_2], times_2)
    events_2[ind] = events_2[ind] + survs[[i]]$n.event[inds_2]
    censors_2[ind] = censors_2[ind] + survs[[i]]$n.censor[inds_2]
  }
  
  risk_1 = rev(rev(c(n[1], n[1]-cumsum(events_1 + censors_1)))[-1])
  risk_2 = rev(rev(c(n[2], n[2]-cumsum(events_2 + censors_2)))[-1])
  
  surv_1 = cumprod(1-(events_1)/risk_1)
  surv_2 = cumprod(1-(events_2)/risk_2)
  
  avgsurv$n = n
  avgsurv$time = c(times_1, times_2)
  avgsurv$n.event = c(events_1, events_2)
  avgsurv$n.censor = c(censors_1, censors_2)
  avgsurv$n.risk = c(risk_1, risk_2)
  avgsurv$surv = c(surv_1, surv_2)
  
  avgsurv$upper = rep(1, length(avgsurv$time))
  avgsurv$lower = rep(0, length(avgsurv$time))
  
  avgsurv$strata[[1]] = length(times_1)
  avgsurv$strata[[2]] = length(times_2)
  
  return (avgsurv)
}


crossvalidatedAUCPlot <- function(predNames, targetname, smoothAUC=FALSE, withNCCN=FALSE, withPSABX=TRUE, withPSA=TRUE, plotSurvivalRatherThanROC = TRUE) {
  
  result <- doCV(X, predNames, targetname, rseed = 1, k_fold=k, ncores=1, smoothAUC=smoothAUC, doPlot=TRUE, plotSurvivalRatherThanROC = plotSurvivalRatherThanROC, col="lightgray", identity.col="darkgray", mark="'", cex.lab=2, cex.axis=1.5)
  if (plotSurvivalRatherThanROC) {
    survs = lapply(result, function(x) x$fit_RF_PSA_SURV)
    lines(avg.surv(survs), mark="'", lwd=5)
    pmodel = median(unlist(lapply(result, function(x) {summary(x$fit_RF_PSA_COX)$logtest[3]})))
  } else {
    rocs = lapply(result, function(x) x$fit_RF_PSA_ROC)
    plot(avg.roc(rocs), add=TRUE, lwd=5)
    cis = mycimedian(unlist(lapply(result, function(x) x$fit_RF_error)))[c(2,3)]
    tmp = unlist(lapply(result, function(x) x$fit_RF_PSA_error))
    text(x=0, y=0.4, labels=paste(c("AUC ", format(median(tmp), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 2, col="black", cex=1.2)
  }
  
  # References PSA
  if (withPSA) {
    resultPSA <- doCV(X, "PSA_Dens", targetname, rseed = 1, k_fold=k, ncores=1, smoothAUC=smoothAUC, doPlot=FALSE, addToPlot = TRUE, col = "dodgerblue4", lty=2)
    if (plotSurvivalRatherThanROC) {
      if (TRUE) {
        survsPSA = lapply(resultPSA, function(x) x$fit_RF_PSA_SURV)
        lines(avg.surv(survsPSA), mark="'", col = "dodgerblue4", lty=2, lwd=5)
        pPSA = median(unlist(lapply(resultPSA, function(x) {summary(x$fit_RF_PSA_COX)$logtest[3]})))
      }
    } else {
      rocsPSA = lapply(resultPSA, function(x) x$fit_RF_PSA_ROC)
      plot(avg.roc(rocsPSA), col = "dodgerblue4", lty=2, lwd=5, add=TRUE)
      cis = mycimedian(unlist(lapply(resultPSA, function(x) x$fit_RF_error)))[c(2,3)]
      tmpPSA = unlist(lapply(resultPSA, function(x) x$fit_RF_PSA_error))
      text(x=0, y=0.3, labels=paste(c("PSA AUC ", format(median(tmpPSA), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 2, col="dodgerblue4", cex=1.2)
    }
  }
  
  # References PSA BX
  if (withPSABX) {
    resultPSABX <- doCV(X, c("PSA_Dens", "BX_GS"), targetname, rseed = 1, k_fold=k, ncores=1, smoothAUC=smoothAUC, doPlot=FALSE, addToPlot = TRUE, col = "firebrick1", lty=4)
    if (plotSurvivalRatherThanROC) {
      if (TRUE) {
        survsPSABX = lapply(resultPSABX, function(x) x$fit_RF_PSA_SURV)
        lines(avg.surv(survsPSABX), mark="'", col = "firebrick1", lty=4, lwd=5)
        pPSABX = median(unlist(lapply(resultPSABX, function(x) {summary(x$fit_RF_PSA_COX)$logtest[3]})))
      }
    } else {
      rocsPSABX = lapply(resultPSABX, function(x) x$fit_RF_PSA_ROC)
      plot(avg.roc(rocsPSABX), col = "firebrick1", lty=4, lwd=5, add=TRUE)
      cis = mycimedian(unlist(lapply(resultPSABX, function(x) x$fit_RF_error)))[c(2,3)]
      tmpPSABX = unlist(lapply(resultPSABX, function(x) x$fit_RF_PSA_error))
      text(x=0, y=0.2, labels=paste(c("PSA/BX AUC ", format(median(tmpPSABX), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 2, col="firebrick1", cex=1.2)
    }
  }
  
  if (withNCCN) {
    if (plotSurvivalRatherThanROC) {
      PredNCCN.surv <- survfit(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ NCCN_risk_group_Bx>2, data = X)
      PredNCCN.cox <-   coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ NCCN_risk_group_Bx>2, data = X)
      lines(PredNCCN.surv, mark="'", lwd=5, lty=5, col="palegreen4")
      pNCCN = summary(PredNCCN.cox)$logtest[3]
    } else {
      rocNCCN = roc(X[, targetname], X$NCCN_risk_group_Bx, smooth = smoothAUC, ci=TRUE)
      plot(rocNCCN, col = "palegreen4", lty=3, lwd=5, add=TRUE)
      cis = rocNCCN$ci[c(1,3)]
      text(x=0, y=0.1, labels=paste(c("NCCN AUC ", format(auc(rocNCCN), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 2, col="palegreen4", cex=1.2)
    }
  }
  
  if (plotSurvivalRatherThanROC) {
    leg = paste0("Model (p=", format(pmodel, digits=3, nsmall=2), ")") 
    leglty = c(1)
    legcol = c("black")
    
    if (withPSA) {
      leg = c(leg, paste0("PSA (p=", format(pPSA, digits=3, nsmall=2), ")"))
      leglty = c(leglty, 2)
      legcol = c(legcol, "dodgerblue4")
    }
    if (withPSABX) {
      leg = c(leg, paste0("PSA/Bx (p=", format(pPSABX, digits=3, nsmall=2), ")"))
      leglty = c(leglty, 4)
      legcol = c(legcol, "firebrick1")
    }
    if (withNCCN) {
      leg = c(leg, paste0("NCCN (p=", format(pNCCN, digits=3, nsmall=2), ")"))
      leglty = c(leglty, 5)
      legcol = c(legcol, "palegreen4")
    }
    legend("bottomleft", leg, col=legcol, lty=leglty, bty="n", cex=2, lwd=5)
  }
  
  return(result)
}


validateSignature = function(predNames, targetname) {

  pdf(paste(c("validate-", paste(predNames, collapse="-"), ".pdf"), collapse=""))
  
  opar = par(no.readonly = TRUE)
  
  result = crossvalidatedAUCPlot(predNames, targetname, plotSurvivalRatherThanROC = FALSE)
  tmp = unlist(lapply(result, function(x) x$fit_RF_PSA_error))
  
  boxes <- cbind(RF_refs_PSA_errors, RF_refs_PSABX_errors, tmp)
  labels <- c("PSA", "PSA + BX", paste(c("PSA", predNames[-1]), collapse=" + "))
  maintit <- "Model compared to PSA and BX"
  par(opar)
  comparisonBoxPlots(lapply(apply(boxes, 2, list), unlist), labels, -1, maintit, ylab="AUC", cex.lab=1.6, cex.axis=1.5, pmar=c(12, 5, 4, 4), cexNames=1)
  
  #for p.values: they have to be adjusted alltogether for the number of tests
  pvals_PSA_adj = p.adjust(c(t.test(tmp, RF_refs_PSA_errors, paired=TRUE)$p.value, pvals_PSA), "BH")
  pval1 = pvals_PSA_adj[1]
  colsub1 = "black"
  if (pval1<=0.05) {
    colsub1 = "red"
  }
  
  pvals_PSABX_adj = p.adjust(c(t.test(tmp, RF_refs_PSABX_errors, paired=TRUE)$p.value, pvals_PSABX), "BH")
  pval2 = pvals_PSABX_adj[1]
  colsub2 = "black"
  if (pval2<=0.05) {
    colsub2 = "red"
  }
  
  mtext(text=paste(c("p=", format(pval1, digits=3, nsmall=2)), collapse=""), col=colsub1, side = 3, at=1, line = 0.2, cex=1)
  mtext(text=paste(c("p=", format(pval2, digits=3, nsmall=2)), collapse=""), col=colsub2, side = 3, at=2, line = 0.2, cex=1)
  
  par(opar)
  
  result = crossvalidatedAUCPlot(predNames, targetname, withPSA = TRUE, plotSurvivalRatherThanROC = TRUE)

  par(opar)
  
  dev.off()

}





for (i in 1:length(predNamesSets)) {
  
  predNames = unique(c("PSA_Dens", sort(predNamesSets[[i]])))
  
  validateSignature(predNames, targetname)
}

  