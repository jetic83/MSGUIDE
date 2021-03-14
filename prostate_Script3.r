source("helpers.r")

delete_successful_lsf(".")





# Directories
wd <- getwd()
td <- sub("/peschuef/", "/peschuef/work/", wd)

### Read all tables
if (tablesNew) {
  if (length(commandArgs()) > 4) {
    z <- as.numeric(commandArgs()[5])
  } else {
    if (!exists("z")) z <- 1
  }
  if (length(commandArgs()) > 5) {
    stepsize <- as.numeric(commandArgs()[6])
  } else {
    if (!exists("stepsize")) stepsize <- 457
  }
  
  GS_RF_PSABX_CVAllz <- vector('list', z)
  
  GS_RF_PSA_CVAllz <- vector('list', z)
  
  GS_RF_CVAllz <- vector('list', z)
  
  PVALS_CVAllz <- vector('list', z)
  
  numberModels <- length(allInd)
  
  for (j in 1:z) {
    GS_RF_PSABX_CVAll <- vector('list', numberModels)
    
    GS_RF_PSA_CVAll <- vector('list', numberModels)
    
    GS_RF_CVAll <- vector('list', numberModels)
    
    PVALS_CVAll <- vector('list', numberModels)
    
    for (i in 1:ceiling(numberModels/stepsize)) {
      tablename <- paste(c(td, "/table", j, "-", i, ".csv"), collapse="")
      show(tablename)
      
      GS_RF_PSABX_CV <- as.list(read.csv2(file=sub("/table", "/RFCV_PSABX", tablename), dec="."))
      
      GS_RF_PSA_CV <- as.list(read.csv2(file=sub("/table", "/RFCV_PSA", tablename), dec="."))
      
      GS_RF_CV <- as.list(read.csv2(file=sub("/table", "/RFCV", tablename), dec="."))
      
      PVALS_CV <- as.list(read.csv2(file=sub("/table", "/PVALS", tablename), dec="."))
      
      
      
      GS_RF_PSABX_CVAll[(((i-1)*stepsize)+1):min(numberModels, (i*stepsize))] <- GS_RF_PSABX_CV
      
      GS_RF_PSA_CVAll[(((i-1)*stepsize)+1):min(numberModels, (i*stepsize))] <- GS_RF_PSA_CV
      
      GS_RF_CVAll[(((i-1)*stepsize)+1):min(numberModels, (i*stepsize))] <- GS_RF_CV
      
      PVALS_CVAll[(((i-1)*stepsize)+1):min(numberModels, (i*stepsize))] <- PVALS_CV
      
      
    }
    GS_RF_PSABX_CVAllz[[j]] <- GS_RF_PSABX_CVAll
    
    GS_RF_PSA_CVAllz[[j]] <- GS_RF_PSA_CVAll
    
    GS_RF_CVAllz[[j]] <- GS_RF_CVAll
    
    PVALS_CVAllz[[j]] <- PVALS_CVAll
    
    
  }
  
  
  
  GS_RF_PSABX_CVAllz <- matrix(GS_RF_PSABX_CVAllz[[1]], nrow=numberModels, ncol=z)
  
  GS_RF_PSA_CVAllz <- matrix(GS_RF_PSA_CVAllz[[1]], nrow=numberModels, ncol=z)
  
  GS_RF_CVAllz <- matrix(GS_RF_CVAllz[[1]], nrow=numberModels, ncol=z)
  
  PVALS_CVAllz <- matrix(PVALS_CVAllz[[1]], nrow=numberModels, ncol=z)
  
  
  # calculate the medians of the errors
  merrorsRF_PSABX <- unlist(lapply(GS_RF_PSABX_CVAllz, function(x){median(x, na.rm=TRUE)}))
  
  merrorsRF_PSA <- unlist(lapply(GS_RF_PSA_CVAllz, function(x){median(x, na.rm=TRUE)}))
  
  merrorsRF <- unlist(lapply(GS_RF_CVAllz, function(x){median(x, na.rm=TRUE)}))
  
  merrorsPVAL <- unlist(lapply(PVALS_CVAllz, function(x){median(x, na.rm=TRUE)}))
  serrorsPVAL <- unlist(lapply(PVALS_CVAllz, function(x){sum(x, na.rm=TRUE)}))
  
  # Order the errors according to their median
  ords <- 1:length(merrorsRF_PSA)
  
  GS_RF_PSABX_CVAll <- GS_RF_PSABX_CVAllz[ords,1]
  
  GS_RF_PSA_CVAll <- GS_RF_PSA_CVAllz[ords,1]
  
  GS_RF_CVAll <- GS_RF_CVAllz[ords,1]
  
  PVALS_CVAll <- PVALS_CVAllz[ords,1]
  
  bsindAll  <- bsindAll[ords]
  merrorsRF_PSABX <- merrorsRF_PSABX[ords]
  
  merrorsRF_PSA <- merrorsRF_PSA[ords]
  
  merrorsRF <- merrorsRF[ords]
  
  merrorsPVAL <- merrorsPVAL[ords]
  serrorsPVAL <- serrorsPVAL[ords]
  
  
  numberModels_notNA <- length(merrorsRF_PSA)
  ind_sig <- 2
  try ({
    while (ind_sig<numberModels_notNA && !all(is.na(GS_RF_PSA_CVAll[[ind_sig]])) && numberModels_notNA*t.test(GS_RF_PSA_CVAll[[ind_sig]], GS_RF_PSA_CVAll[[1]], paired=TRUE)$p.value >= 0.05) {
      ind_sig <- ind_sig+1
    }
  } , silent=TRUE)
  
  
  tablesNew <- FALSE
  
  save.image()
  save.image(file = paste(c(gsub("/peschuef/", "/peschuef/work/", getwd()), "/", pdfname, ".RData"), collapse = ""))
  
}

# PDF
pdffilename <- paste(c(pdfname, "_", valtyp, "-", as.character(k), "_", as.character(seed_script1)), collapse="")
pdf(paste(c(pdffilename, ".pdf"), collapse=""))
par_orig <- par(no.readonly = TRUE)


# Find a specific model
ind = findModelIndex(c("VTNC", "FINC"), predNames20, bsindAll)







###
### Helper Function plot comparison boxes
###
# PLOT THE BOXES
comparisonBoxPlots <- function(boxes, labels, n_models, maintit, cex.axis=0.7, pmar = c(12, 4, 4, 4), cexNames=0.7, ...) {
  par(mar=pmar)
  myboxplot(boxes,
            col=c("lightblue", "blue", rep("orange", length(labels)-2)), xaxt = "n",  xlab = "",
            las=2, cex.axis=cex.axis, main=maintit, notch=FALSE, varwidth=FALSE, ...)
  # x axis with ticks but without labels
  axis(1, labels = FALSE, at=seq_along(labels))
  # Plot x labs at default x position
  text(x =  seq_along(labels), y = par("usr")[3] - 0.02, srt = 40, adj = 1,
       labels = labels, xpd = TRUE, cex=cexNames)
  if (!is.na(n_models) && n_models>=0) {
    mtext(paste(as.character(n_models), "models", sep=" "), 1, line=10)
  }
  # end PLOT THE BOXES
}

###
### Helper Function plot comparison boxes
###
# PLOT THE BOXES
comparisonBoxPlotsWOReferences <- function(boxes, labels, n_models, maintit) {
  par(mar=c(12, 4, 4, 4))
  myboxplot(boxes,
            col="orange", xaxt = "n",  xlab = "",
            las=2, cex.axis=0.7, main=maintit, notch=FALSE, varwidth=FALSE)
  # x axis with ticks but without labels
  axis(1, labels = FALSE, at=seq_along(labels))
  # Plot x labs at default x position
  text(x =  seq_along(labels), y = par("usr")[3] - 0.02, srt = 40, adj = 1,
       labels = labels, xpd = TRUE, cex=0.7)
  mtext(paste(as.character(n_models), "models", sep=" "), 1, line=10)
  # end PLOT THE BOXES
}

###
### Helper Function plot comparison AUC
###
# PLOT THE AUC
comparisonAUCPlotsWOReferences <- function(predNames, dat_train, dat_test, r_seed = rseed, ctype="rf", ...) {
  set.seed(r_seed)
  if(ctype=="rf") {
    fit_RF <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c(predNames)), collapse="+")), collapse="")), data.frame(dat_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(unique(c(predNames))))
    prediction = predict(fit_RF, newdata = dat_test, type="prob")[,2]
  } else {
    fit_RF <-          glm(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c(predNames)), collapse="+")), collapse="")), family=binomial, data=data.frame(dat_train))
    prediction = predict(fit_RF, newdata = dat_test, type="response")
  }
  
  roc_ = roc(dat_test[, targetname], prediction)
  cutoff = roc_$thresholds[which.max(roc_$specificities + roc_$sensitivities)[1]]
  
  return (c(errormeasureAUC(dat_test[, targetname], prediction, doPlot = TRUE, ...),
            cutoff))
}

comparisonAUCPlotsWReferences <- function(predNames, dat_train, dat_test, r_seed = rseed, ctype="rf", ...) {
  par(par_orig)
  res <- comparisonAUCPlotsWOReferences(predNames, dat_train, dat_test, r_seed = rseed, add=FALSE, print.thres=TRUE, ctype=ctype, ...)
  auc1 <- res[1]
  cutoff <- res[2] 
  res <- comparisonAUCPlotsWOReferences(c("PSA_Dens"), dat_train, dat_test, r_seed = rseed, col="red", lty=2, add=TRUE, ctype=ctype, ...)
  auc2 <- res[1]
  cutoffPSA <- res[2] 
  res <- comparisonAUCPlotsWOReferences(c("PSA_Dens", "BX_GS"), dat_train, dat_test, r_seed = rseed, col="orange", lty=3, add=TRUE, ctype=ctype, ...)
  auc3 <- res[1]
  cutoffPSABX <- res[2] 
  
  roc_ = roc(dat_test$BCR5y, dat_test$NCCN_risk_group_Bx)
  auc4 <- auc(roc_)
  plot(roc_, add=TRUE, col="green", lty=4)
  
  title(main=paste(predNames, collapse=" + "))
  legend(x = 0.4, y=0.2, legend=c(paste("AUC", format(auc1, digits=3, nsmall=2)), paste("PSA AUC", format(auc2, digits=3, nsmall=2)), paste("PSA/BX AUC", format(auc3, digits=3, nsmall=2)), paste("NCCN AUC", format(auc4, digits=3, nsmall=2))), lty=1:4, col=c("black", "red", "orange", "green"), bty="n")
  
  return (c(cutoff, cutoffPSA, cutoffPSABX))
}

# Median error of reference Models PSA and PSABX
merror_RF_refs_PSA <- median(RF_refs_PSA_errors)
merror_RF_refs_PSABX <- median(RF_refs_PSABX_errors)

# Find best RF models better than references
ind_models_RF_PSA_best <- which(merrorsRF_PSA >= merror_RF_refs_PSABX & merrorsRF_PSA >= merror_RF_refs_PSA)

# Record number of best RF models
n_RF_PSA_best <- length(ind_models_RF_PSA_best)

# sort best RF models 
ind = order(merrorsRF_PSA[ind_models_RF_PSA_best], decreasing = TRUE)
ind_models_RF_PSA_best = ind_models_RF_PSA_best[ind]

# Find the first models not differing significantly.
ind_RF_PSA_sig <- 2
try ({
  while (ind_RF_PSA_sig < numberModels_notNA && !all(is.na(GS_RF_PSA_CVAll[[ind_models_RF_PSA_best[ind_RF_PSA_sig]]])) && numberModels_notNA*t.test(GS_RF_PSA_CVAll[[ind_models_RF_PSA_best[ind_RF_PSA_sig]]], GS_RF_PSA_CVAll[[ind_models_RF_PSA_best[1]]], paired=TRUE)$p.value >= 0.05) {
    ind_RF_PSA_sig <- ind_RF_PSA_sig+1
  }
} , silent=TRUE)


if (n_RF_PSA_best>0) {
  predHist2(predNames20[unlist(bsindAll[ind_models_RF_PSA_best[1:ind_RF_PSA_sig]])], tit="First non differing models better than PSA+BX")
  
  ## Rankings over Top Models
  nmin <- 1
  nmax <- 2*ind_RF_PSA_sig
  max_res <- 500
  sampling <- if (nmax-nmin+1<max_res) nmin:nmax else round(seq(nmin, nmax, (nmax-nmin+1)/max_res))
  sampling[length(sampling)] <- nmax
  H <- t(sapply(sampling, function(y){
    distr <- rep(0, length(predNames20))
    names(distr) = predNames20
    ph = predHist2(predNames20[unlist(bsindAll[ind_models_RF_PSA_best[1:y]])], doPlot=FALSE, normalized=FALSE, n=y)
    distr[names(ph)] = ph
    distr
  }))
  par(mai=c(1.02, 0.86, 0.82, 0.02))
  plot(c(0,1.4*nmax), c(0, 1), type="n", xlab="Top n Models", ylab="Frequency of Features [%]", main="Ranking of Features in Top Models", cex.main=1.7, cex.axis=1.5, cex.lab=1.5)
  if (ind_RF_PSA_sig<=nmax) abline(v=ind_RF_PSA_sig, col="grey")
  abline(h=0.5, col="lightgrey")
  for (i in 1:dim(H)[2]) lines(sampling, H[,i], col=rainbow(dim(H)[2])[i])
  text(nmax, H[dim(H)[1],], labels=colnames(H), col=rainbow(dim(H)[2]), pos=4, cex=1)
  par(par_orig)
  ## END Rankings over Top Models
  

  par(mar=c(12, 4, 4, 4))
  
  # PLOT THE BOXES
  boxes <- cbind(RF_refs_PSA_errors, RF_refs_PSABX_errors, sapply(ind_models_RF_PSA_best[1:min(10,n_RF_PSA_best)], function(x) GS_RF_PSA_CVAll[[x]]))
  labels <- c("PSA", "PSA + BX", lapply(lapply(ind_models_RF_PSA_best[1:min(10,n_RF_PSA_best)], function(x) predNames20[unlist(bsindAll[x])]), function(y) paste(y, collapse=" + ")))
  maintit <- "Best RF (incl PSA) Models compared to PSA and BX"
  comparisonBoxPlots(lapply(apply(boxes, 2, list), unlist), labels, n_RF_PSA_best, maintit)
}



# find the models which are better than the references per fold
ref_max <- sapply(1:k, function(y) max(RF_refs_PSABX_errors[y])) 

GS_RF_CVAll_pf <- lapply(GS_RF_CVAll, function(x) (x - ref_max))

merrorsRF_pf <- unlist(lapply(GS_RF_CVAll_pf, function(x){median(x, na.rm=TRUE)}))

# Find best RF models better than references in each fold
ind_models_RF_best_pf <- which(merrorsRF_pf > 0)

# Record number of best RF models (best per fold)
n_RF_best_pf <- length(ind_models_RF_best_pf)

# sort best RF models 
ind = order(merrorsRF_pf[ind_models_RF_best_pf], decreasing = TRUE)
ind_models_RF_best_pf = ind_models_RF_best_pf[ind]


# Find the first models not differing significantly.
ind_RF_sig_pf <- 2
try ({
  while (ind_RF_sig_pf <= length(ind_models_RF_best_pf) && !all(is.na(GS_RF_CVAll_pf[[ind_models_RF_best_pf[ind_RF_sig_pf]]])) && numberModels_notNA*t.test(GS_RF_CVAll_pf[[ind_models_RF_best_pf[ind_RF_sig_pf]]], GS_RF_CVAll_pf[[ind_models_RF_best_pf[1]]], paired=TRUE)$p.value >= 0.05) {
    ind_RF_sig_pf <- ind_RF_sig_pf+1
  }
} , silent=TRUE)


if (n_RF_best_pf>0) {
  predHist2(predNames20[unlist(bsindAll[ind_models_RF_best_pf[1:ind_RF_sig_pf]])], tit="First non differing models better than PSA+BX")
  
  ## Rankings over Top Models
  nmin <- 1
  nmax <- 2*ind_RF_sig_pf
  max_res <- 500
  sampling <- if (nmax-nmin+1<max_res) nmin:nmax else round(seq(nmin, nmax, (nmax-nmin+1)/max_res))
  sampling[length(sampling)] <- nmax
  H <- t(sapply(sampling, function(y){
    distr <- rep(0, length(predNames20))
    names(distr) = predNames20
    ph = predHist2(predNames20[unlist(bsindAll[ind_models_RF_best_pf[1:y]])], doPlot=FALSE, normalized=FALSE, n=y)
    distr[names(ph)] = ph
    distr
  }))
  par(mai=c(1.02, 0.86, 0.82, 0.02))
  plot(c(0,1.4*nmax), c(0, 1), type="n", xlab="Top n Models", ylab="Frequency of Features [%]", main="Ranking of Features in Top Models", cex.main=1.7, cex.axis=1.5, cex.lab=1.5)
  if (ind_RF_sig_pf<=nmax) abline(v=ind_RF_sig_pf, col="grey")
  abline(h=0.5, col="lightgrey")
  for (i in 1:dim(H)[2]) lines(sampling, H[,i], col=rainbow(dim(H)[2])[i])
  text(nmax, H[dim(H)[1],], labels=colnames(H), col=rainbow(dim(H)[2]), pos=4, cex=1)
  par(par_orig)
  ## END Rankings over Top Models

  
  # PLOT THE BOXES
  boxes <- cbind(sapply(ind_models_RF_best_pf[1:min(10,n_RF_best_pf)], function(x) GS_RF_CVAll_pf[[x]]))
  labels <- c(lapply(lapply(ind_models_RF_best_pf[1:min(10,n_RF_best_pf)], function(x) predNames20[unlist(bsindAll[x])]), function(y) paste(y, collapse=" + ")))
  maintit <- "Best RF Models compared to PSA and BX per fold"
  comparisonBoxPlotsWOReferences(lapply(apply(boxes, 2, list), unlist), labels, n_RF_best_pf, maintit)
  
}


dev.off()


save.image(file = paste(c(gsub("/peschuef/", "/peschuef/work/", getwd()), "/", pdfname, ".RData"), collapse = ""))


pvals_PSA = sapply(1:numberModels_notNA, function(x) t.test(GS_RF_PSA_CVAll[[x]], RF_refs_PSA_errors, paired=TRUE)$p.value)
pvals_PSABX = sapply(1:numberModels_notNA, function(x) t.test(GS_RF_PSABX_CVAll[[x]], RF_refs_PSA_errors, paired=TRUE)$p.value)
ind = findModelIndex(c("VTNC", "FINC"), predNames20, bsindAll)
p.adjust(pvals_PSA, "BH")[ind]
p.adjust(pvals_PSABX, "BH")[ind]

source("validation.r", verbose=TRUE)

