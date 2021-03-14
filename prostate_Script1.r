
source("helpers.r")


# Parameter setting
if (length(commandArgs()) > 4) {
  seed_script1 <- as.numeric(commandArgs()[5])
} else {
  if (!exists("seed_script1")) seed_script1 <- 1
}

if (length(commandArgs()) > 5) {
  batchj <- as.numeric(commandArgs()[6])
} else {
  if (!exists("batchj")) batchj <- 4
}



rseed = seed_script1

####
####
#### Read the Data
####
alldata = read.csv("alldata.csv", sep = ",")

source("helpers.r")

pdfname <- "UNSC" # unscaled
targetname <- "BCR5y"
valtyp <- "LOOCV" 
errorfun <- "aucerror"
cormethod <- "spearman" 
cl_family <- "binomial" 
k <- 50  # number of bootsraps per run and model. k=100 means 100 BSs, so 100 errors, formed to one boxplot with one median
z <- 1    # number of runs. z=30 means 30 runs, so histograms are averaged over 30 histograms, there are 30 distances of histograms to the average histogram.
e <- 1000 # end. Look at the e best models at most (default 1000).
ntree <- 200 # number of trees for a random forest
ctype="rf"



predNames20 = c("ATRN", "BTD", "CADM1", "CATD", "CERU", "CFAH", "ECM1", "FINC", "GOLM1", "HYOU1", "ICAM1", "LG3BP",
                "LUM", "NCAM1", "PIGR", "PLXB2", "POSTN", "TRFE", "TSP1", "VTNC", "ZA2G")

predNames20 <- sort(predNames20)

predNamesRef_PSA <- c("PSA_Dens")
predNamesRef_PSABX <- c(predNamesRef_PSA, "BX_GS")


alldata <- removeNas(alldata, unique(c(predNamesRef_PSABX, predNames20, targetname)))




# Train on one batch ProCOC X
train_data_global <- alldata[which(alldata$reading==1 & alldata$ProCOC==batchj), ]
test_data_global <- alldata[which(alldata$reading==1 & alldata$ProCOC!=batchj), ]
pdfname <- paste(c(pdfname, "_batch", batchj), collapse = "")
pdfname <- paste(c(pdfname, "_", errorfun, "_", ctype), collapse="")

dat_train <- train_data_global
dat_test  <- test_data_global # for overall validation, will not be used in CV

dat_train <- removeNas(dat_train, unique(c(predNamesRef_PSABX, predNames20, targetname)))
dat_test  <- removeNas(dat_test, unique(c(predNamesRef_PSABX, predNames20, targetname)))

X <- dat_train


n <- length(predNames20)  # number of predictors from which a combination can be chosen (default 20)
p <- 5  # highest number of predictors in one model (default 5)

# LOAD INDICES
indexfile <- paste(as.character(n), as.character(p), ".RData", sep="")
if (file.exists(indexfile)) {
  load(indexfile)
} else {
  allInd <- CreateIndices(n, p) 
  save(allInd, file=indexfile)
}



# BOOTSTRAPS            
indRowNA <- dim(X)[1]+1
d = dim(X[-indRowNA,])[1]
allind <- c(1:d)
bserrorsAllz <- list()
bsindAll <- allInd

tablesNew <- TRUE



# calculate reference models PSA
RF_refs_PSABX <- list()
GLM_refs_PSABX <- list()
RF_refs_PSA <- list()
GLM_refs_PSA <- list()

RF_refs_PSABX_errors <- c()
GLM_refs_PSABX_errors <- c()
RF_refs_PSA_errors <- c()
GLM_refs_PSA_errors <- c()

train_inds <- list()
test_inds <- list()
set.seed(rseed)
for (i in 1:k) {
  set.seed(i)
  ind <- order(dat_train[,targetname])[sampling::strata(dat_train[order(dat_train[,targetname]),], targetname, round(0.69 * as.numeric(summary(dat_train[, targetname])[2:3])), "srswor")$ID_unit]
  ind_oob <- (1:nrow(dat_train))[-ind]
  
  train_inds <- c(train_inds, list(ind))
  test_inds <- c(test_inds, list(ind_oob))
  
  dat_fold_train <- dat_train[ind, ]
  dat_fold_test  <- dat_train[ind_oob, ]
  
  
  ##
  ## Random Forest
  ##
  require(randomForest)
  set.seed(rseed)
  rf_fit <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(predNamesRef_PSABX, collapse="+")), collapse="")), data.frame(dat_fold_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(predNamesRef_PSABX))
  #RF_refs_PSABX_errest <- c(RF_refs_PSABX_errest, 1-sum(diag(rf_fit$confusion))/sum(rf_fit$confusion[1:2, 1:2]))
  if (errorfun == "aucerror") {
    RF_refs_PSABX_errors <- c(RF_refs_PSABX_errors, errormeasureAUC(dat_fold_test[, targetname], predict(rf_fit, newdata = dat_fold_test, type="prob")[,2]))
  } else if (errorfun == "surverror") {
    RF_refs_PSABX_errors <- c(RF_refs_PSABX_errors, errormeasureSURV(predict(rf_fit, newdata = dat_fold_test), dat_fold_test))
  } else {
    RF_refs_PSABX_errors <- c(RF_refs_PSABX_errors, errormeasure(dat_fold_test[, targetname], predict(rf_fit, newdata = dat_fold_test)))
  }
  RF_refs_PSABX <- c(RF_refs_PSABX, list(rf_fit))
  
  set.seed(rseed)
  rf_fit <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(predNamesRef_PSA, collapse="+")), collapse="")), data.frame(dat_fold_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(predNamesRef_PSA))
  if (errorfun == "aucerror") {
    RF_refs_PSA_errors <- c(RF_refs_PSA_errors, errormeasureAUC(dat_fold_test[, targetname], predict(rf_fit, newdata = dat_fold_test, type="prob")[,2]))
  } else if (errorfun == "surverror") {
    RF_refs_PSA_errors <- c(RF_refs_PSA_errors, errormeasureSURV(predict(rf_fit, newdata = dat_fold_test), dat_fold_test))
  } else {
    RF_refs_PSA_errors <- c(RF_refs_PSA_errors, errormeasure(dat_fold_test[, targetname], predict(rf_fit, newdata = dat_fold_test)))
  }
  RF_refs_PSA <- c(RF_refs_PSA, list(rf_fit))
  
}



###
### Validation Functions
###

doClassifierSurvivalFitWithClassifier <- function(data, classifier, tit, classifierPSA=NA, classifierPSABX=NA, legend=TRUE, doPlot=TRUE, plotSurvivalRatherThanROC=TRUE, plotPSA=TRUE, plotPSABX=TRUE, plotNCCN=TRUE, cutoff=0.5, cutoffPSA=0.5, cutoffPSABX=0.5, ctype="rf", mark = 3, lwd=1, ...) {
  require(pROC)
  require(survival)
  
  data. <- data
  
  if (ctype=="rf") {
    data.$Score <- predict(classifier, newdata = data., type="prob")[,2]
  } else {
    data.$Score <- predict(classifier, newdata = data., type="response")
  }
  
  rocScore = roc(data.[, targetname], data.$Score, ci=TRUE)
  if (cutoff<0) {
    cutoff = rocScore$thresholds[which.max(rocScore$specificities + rocScore$sensitivities)[1]]
  }
  if (class(data.$Score) == "numeric") {
    data.$Score <- factor(data.$Score >= cutoff)
  }
  levels(data.$Score) <- c("LOW", "HIGH")
  
  Pred.surv <- survfit(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ Score, data = data.)
  Pred.cox  <-   coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ Score, data = data.)
  
  if (doPlot) {
    if (plotSurvivalRatherThanROC)
      plot(Pred.surv, col = 1:2, ylab="Cumulative Survival", xlab="Months", main=tit, sub=paste(c("n = ", paste(Pred.surv$n, collapse=" / "), ",  p=", format(summary(Pred.cox)$logtest[3], digits=3, nsmall=2)), collapse=""), lwd=lwd, mark=mark, ...)
    else {
      plot(rocScore, col = 1, main=tit, sub=paste(c("n = ", length(rocScore$cases)+length(rocScore$controls)), collapse=""), mar=c(5.1,4.1,4.2,2.1), ...)
      cis = rocScore$ci[c(1,3)]
      text(x=0.4, y=0.4, labels=paste(c("AUC ", format(auc(rocScore), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 4, col=1)
    }
    
    if (plotPSA && !is.na(classifierPSA)) {
      if (ctype=="rf") {
        data.$PSAScore <- predict(classifierPSA, newdata = data., type="prob")[,2]
      } else {
        data.$PSAScore <- predict(classifierPSA, newdata = data., type="response")
      }
      rocPSA = roc(data.[, targetname], data.$PSAScore, ci=TRUE)
      if (cutoffPSA<0) {
        cutoffPSA = rocPSA$thresholds[which.max(rocPSA$specificities + rocPSA$sensitivities)[1]]
      }
      if (class(data.$PSAScore) == "numeric") {
        data.$PSAScore <- factor(data.$PSAScore >= cutoffPSA)
      }
      levels(data.$PSAScore) <- c("LOW", "HIGH")
      
      PredPSA.surv <- survfit(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ PSAScore, data = data.)
      PredPSA.cox  <-   coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ PSAScore, data = data.)
      
      if (plotSurvivalRatherThanROC) {
        lines(PredPSA.surv, col = c("gray", "orange"), lty=2, mark=mark, lwd=lwd)
        mtext(text=paste(c("PSA ", paste(PredPSA.surv$n, collapse="/"), " (p=", format(summary(PredPSA.cox)$logtest[3], digits=3, nsmall=2), ")"), collapse=""), col="darkgray", side = 3, adj = 0, cex=0.7)
      } else {
        plot(rocPSA, col = "darkgray", lty=2, add=TRUE, ...)
        mtext(text=paste(c("PSA ", paste(c("n = ", length(rocPSA$cases)+length(rocPSA$controls)), collapse=""), " (p=", format(roc.test(rocScore, rocPSA)$p.value, digits=3, nsmall=2), ")"), collapse=""), col="darkgray", side = 3, adj = 0, cex=0.7)
        cis = rocPSA$ci[c(1,3)]
        text(x=0.4, y=0.3, labels=paste(c("AUC ", format(auc(rocPSA), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 4, col="darkgray")
      }
    }
    
    if (plotPSABX && !is.na(classifierPSABX)) {
      if (ctype=="rf") {
        data.$PSABXScore <- predict(classifierPSABX, newdata = data., type="prob")[,2]
      } else {
        data.$PSABXScore <- predict(classifierPSABX, newdata = data., type="response")
      }
      rocPSABX = roc(data.[, targetname], data.$PSABXScore, ci=TRUE)
      if (cutoffPSABX<0) {
        cutoffPSABX = rocPSABX$thresholds[which.max(rocPSABX$specificities + rocPSABX$sensitivities)[1]]
      }
      if (class(data.$PSABXScore) == "numeric") {
        data.$PSABXScore <- factor(data.$PSABXScore >= cutoffPSABX)
      }
      levels(data.$PSABXScore) <- c("LOW", "HIGH")
      
      PredPSABX.surv <- survfit(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ PSABXScore, data = data.)
      PredPSABX.cox  <-   coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ PSABXScore, data = data.)
      
      if (plotSurvivalRatherThanROC) {
        lines(PredPSABX.surv, col = c("lightblue3", "lightblue4"), lty=4, mark=mark, lwd=lwd)
        mtext(text=paste(c("PSA/BX ", paste(PredPSABX.surv$n, collapse="/"), " (p=", format(summary(PredPSABX.cox)$logtest[3], digits=3, nsmall=2), ")"), collapse=""), col="lightblue3", side = 3, adj = 1, cex=0.7)
      } else {
        plot(rocPSABX, col = "lightblue3", lty=4, add=TRUE, ...)
        mtext(text=paste(c("PSA/BX ", paste(c("n = ", length(rocPSABX$cases)+length(rocPSABX$controls)), collapse=""), " (p=", format(roc.test(rocScore, rocPSABX)$p.value, digits=3, nsmall=2), ")"), collapse=""), col="lightblue3", side = 3, adj = 1, cex=0.7)
        cis = rocPSABX$ci[c(1,3)]
        text(x=0.4, y=0.2, labels=paste(c("AUC ", format(auc(rocPSABX), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 4, col="lightblue3")
      }
    }
    
    if (plotNCCN) {
      # References NCCN
      rocNCCN = roc(data.[, targetname], data.$NCCN_risk_group_Bx, ci=TRUE)
      PredNCCN.surv <- survfit(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ NCCN_risk_group_Bx>2, data = data.)
      PredNCCN.cox  <-   coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ NCCN_risk_group_Bx>2, data = data.)
      
      if (plotSurvivalRatherThanROC) {
        lines(PredNCCN.surv, col = c("palegreen3", "palegreen4"), lty=5, mark=mark, lwd=lwd)
        mtext(text=paste(c("NCCN ", paste(PredNCCN.surv$n, collapse="/"), " (p=", format(summary(PredNCCN.cox)$logtest[3], digits=3, nsmall=2), ")"), collapse=""), col="palegreen3", side = 4, adj = 1, cex=0.7)
      } else {
        plot(rocNCCN, col = "palegreen3", lty=5, add=TRUE, ...)
        mtext(text=paste(c("NCCN ", paste(c("n = ", length(rocNCCN$cases)+length(rocNCCN$controls)), collapse=""), " (p=", format(roc.test(rocScore, rocNCCN)$p.value, digits=3, nsmall=2), ")"), collapse=""), col="palegreen3", side = 4, adj = 1, cex=0.7)
        cis = rocNCCN$ci[c(1,3)]
        text(x=0.4, y=0.1, labels=paste(c("AUC ", format(auc(rocNCCN), digits=3, nsmall=2), " [", format(cis[1], digits=3, nsmall=2), ";", format(cis[2], digits=3, nsmall=2), "]"), collapse=""), pos = 4, col="palegreen3")
      }
    }
    
    if (legend) {
      if (plotSurvivalRatherThanROC && !is.null(names(Pred.surv$strata))) {
        leg = names(Pred.surv$strata)
        leglty = c(1,1)
        legcol = c("black", "red")
        
        if (plotPSA && !is.na(classifierPSA)) {
          leg = c(leg, names(PredPSA.surv$strata))
          leglty = c(leglty, c(2,2))
          legcol = c(legcol, c("gray", "orange"))
        }
        if (plotPSABX && !is.na(classifierPSABX)) {
          leg = c(leg, names(PredPSABX.surv$strata))
          leglty = c(leglty, c(4,4))
          legcol = c(legcol, c("lightblue3", "lightblue4"))
        }
        if (plotNCCN) {
          leg = c(leg, names(PredNCCN.surv$strata))
          leglty = c(leglty, c(5,5))
          legcol = c(legcol, c("palegreen3", "palegreen4"))
        }
        legend("bottomleft", leg, col=legcol, lty=leglty, bty="n", cex=0.8)
      }
    }
  }
  
  return (summary(Pred.cox)$logtest[3])
}

doClassifierSurvivalFitWithPrednames <- function(datatrain, datatest, predNames, targetname, tit="", ranseed = rseed, ktree = ntree, doPlot=TRUE, cutoff=0.5, cutoffPSA=0.5, cutoffPSABX=0.5, ctype="rf", plotSurvivalRatherThanROC=TRUE, ...) {
  
  set.seed(ranseed)
  if (ctype=="rf") {
    require(randomForest)
    classifier = randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(unique(predNames), collapse="+")), collapse="")), data.frame(datatrain), importance=TRUE, na.action=na.omit, ntree = ktree, mtry=length(unique(predNames)))
  } else {
    classifier = glm(formula(paste(c("as.factor(", targetname, ")~", paste(unique(predNames), collapse="+")), collapse="")), family=binomial, data=data.frame(datatrain))
  }
  
  set.seed(ranseed)
  if (ctype=="rf") {
    classifierPSA = randomForest(formula(paste(c("as.factor(", targetname, ")~", paste("PSA_Dens", collapse="+")), collapse="")), data.frame(datatrain), importance=TRUE, na.action=na.omit, ntree = ktree)
  } else {
    classifierPSA = glm(formula(paste(c("as.factor(", targetname, ")~", paste("PSA_Dens", collapse="+")), collapse="")), family=binomial, data=data.frame(datatrain))
  }
  
  set.seed(ranseed)
  if (ctype=="rf") {
    classifierPSABX = randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c("PSA_Dens", "BX_GS")), collapse="+")), collapse="")), data.frame(datatrain), importance=TRUE, na.action=na.omit, ntree = ktree, mtry=2)
  } else {
    classifierPSABX = glm(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c("PSA_Dens", "BX_GS")), collapse="+")), collapse="")), family=binomial, data=data.frame(datatrain))
  }
  
  if (is.na(tit) || tit=="") {
    tit = paste(predNames, collapse=" + ")
  }
  
  return (doClassifierSurvivalFitWithClassifier(datatest, classifier, tit, classifierPSA, classifierPSABX, doPlot=doPlot, cutoff=cutoff, cutoffPSA=cutoffPSA, cutoffPSABX=cutoffPSABX, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...))
}

doSurvivalTestWithPrednames <- function(predNames, alldata.=alldata, doPlot=TRUE, ctype="rf", plotSurvivalRatherThanROC=TRUE, ...) {
  par_before <- par(no.readonly = TRUE)
  par(mfrow=c(2,2))
  I   <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC!=1), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==1), ], predNames, targetname, tit="ProCOC I", legend=FALSE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
  II  <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC!=2), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==2), ], predNames, targetname, tit="ProCOC II", legend=FALSE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
  III <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC!=3), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==3), ], predNames, targetname, tit="ProCOC III", legend=FALSE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
  IV  <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC!=4), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==4), ], predNames, targetname, tit="ProCOC IV", legend=TRUE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
  if (doPlot) {
    at. = 1.5
    if (plotSurvivalRatherThanROC)
      at. = -10
    mtext(paste(predNames, collapse=" + "),line = 24, at = at., col="sienna3")
    mtext(paste(c("Batchwise\nCross-Validated"), collapse=""),line = 3.5, at = at., col="tomato")
  }
  par(par_before)
  return (c(I=I, II=II, III=III, IV=IV))
}

doSurvivalTestWithPrednamesBatchwise <- function(predNames, alldata.=alldata, doPlot=TRUE, ctype="rf", plotSurvivalRatherThanROC=TRUE, ...) {
  par_before <- par(no.readonly = TRUE)
  par(mfrow=c(2,2))
  allbatches <- c()
  for (batch in 1:4) {
    I   <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC==batch), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==1), ], predNames, targetname, tit="ProCOC I", legend=FALSE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
    II  <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC==batch), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==2), ], predNames, targetname, tit="ProCOC II", legend=FALSE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
    III <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC==batch), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==3), ], predNames, targetname, tit="ProCOC III", legend=FALSE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
    IV  <- doClassifierSurvivalFitWithPrednames(alldata.[which(alldata.$reading==1 & alldata.$ProCOC==batch), ], alldata.[which(alldata.$reading==1 & alldata.$ProCOC==4), ], predNames, targetname, tit="ProCOC IV", legend=TRUE, doPlot=doPlot, cutoff=-1, cutoffPSA=-1, cutoffPSABX=-1, ctype=ctype, plotSurvivalRatherThanROC=plotSurvivalRatherThanROC, ...)
    if (doPlot) {
      at. = 1.5
      if (plotSurvivalRatherThanROC)
        at. = -10
      mtext(paste(predNames, collapse=" + "),line = 24, at = at., col="sienna3")
      mtext(paste(c("Trained on Batch ", batch), collapse=""),line = 3.5, at = at., col="tomato")
    }
    allbatches <- c(allbatches, list(c(I=I, II=II, III=III, IV=IV)))
  }
  par(par_before)
  return (allbatches)
}

doCV <- function(dat_train, predNames, targetname, rseed = 1, testdata = NA, k_fold = k, ncores = round(0.5*detectCores()), doPlot=FALSE, addToPlot=NA, tit=NA, smoothAUC=FALSE, plotSurvivalRatherThanROC=TRUE, ...) {
  require(doParallel)
  require(pROC)
  require(survival)
  require(randomForestSRC)
  
  registerDoParallel(cores=ncores)
  
  
  CVfolds <- foreach(i=1:k_fold,.combine="c", .packages = c('randomForest', 'boot', 'survival', 'randomForestSRC')) %do% { #inner loop. Change %dopar% to %do% to run sequentially (a normal for loop)
    
    if (is.na(addToPlot)) {
      addToPlot_ = i!=1
    } else {
      addToPlot_ = addToPlot
    }
    if (is.na(tit)) {
      tit_ = paste(predNames, collapse=" + ")
    } else {
      tit_ = tit
    }
    
    ind <- train_inds[[i]]
    ind_oob <- test_inds[[i]]
    
    dat_fold_train <- dat_train[ind, ]
    dat_fold_test  <- dat_train[ind_oob, ]
    
    CVfold <- c()
    
    ##
    ## Random Forest
    ##
    if (ctype=="rf") {
      require(randomForest)
      set.seed(rseed)
      CVfold$fit_RF_PSABX       <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c(predNamesRef_PSABX, predNames)), collapse="+")), collapse="")), data.frame(dat_fold_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(unique(c(predNamesRef_PSABX, predNames))))
      if (errorfun == "aucerror") {
        CVfold$fit_RF_PSABX_error <- errormeasureAUC(dat_fold_test[, targetname], predict(CVfold$fit_RF_PSABX, newdata = dat_fold_test, type="prob")[,2])
      } 
       
      set.seed(rseed)
      CVfold$fit_RF_PSA       <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c(predNamesRef_PSA, predNames)), collapse="+")), collapse="")), data.frame(dat_fold_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(unique(c(predNamesRef_PSA, predNames))))
      if (errorfun == "aucerror") {
        CVfold$fit_RF_PSA_error <- errormeasureAUC(dat_fold_test[, targetname], predict(CVfold$fit_RF_PSA, newdata = dat_fold_test, type="prob")[,2], doPlot=(doPlot && !plotSurvivalRatherThanROC), smoothAUC, add=addToPlot_, main=tit_, getROC = FALSE, ...)
        CVfold$fit_RF_PSA_ROC   <- errormeasureAUC(dat_fold_test[, targetname], predict(CVfold$fit_RF_PSA, newdata = dat_fold_test, type="prob")[,2], doPlot=FALSE, smoothAUC, add=addToPlot_, main=tit_, getROC = TRUE, ...)
        
        if (!smoothAUC) {
          roc_ = CVfold$fit_RF_PSA_ROC
          cutoff = roc_$thresholds[which.max(roc_$specificities + roc_$sensitivities)[1]]
          dat_fold_test$Score = predict(CVfold$fit_RF_PSA, newdata = dat_fold_test, type="prob")[,2]
          if (class(dat_fold_test$Score) == "numeric") {
            dat_fold_test$Score <- factor(dat_fold_test$Score >= cutoff)
          }
          levels(dat_fold_test$Score) <- c("LOW", "HIGH")
          
          Pred.surv <- survfit(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ Score, data = dat_fold_test)
          Pred.cox  <-   coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ Score, data = dat_fold_test)
          
          CVfold$fit_RF_PSA_SURV <- Pred.surv
          CVfold$fit_RF_PSA_COX <- Pred.cox
          
          if (doPlot && plotSurvivalRatherThanROC) {
            if (!addToPlot_) {
              plot(Pred.surv, ylab="Cumulative Survival", xlab="Months", main=tit_, ...)
            } else {
              lines(Pred.surv, ...)
            }
          }
        }
      }
      
      set.seed(rseed)
      CVfold$fit_RF       <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c(predNames)), collapse="+")), collapse="")), data.frame(dat_fold_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(unique(c(predNames))))
      if (errorfun == "aucerror") {
        CVfold$fit_RF_error <- errormeasureAUC(dat_fold_test[, targetname], predict(CVfold$fit_RF, newdata = dat_fold_test, type="prob")[,2], getROC=FALSE)
        CVfold$fit_RF_ROC <- errormeasureAUC(dat_fold_test[, targetname], predict(CVfold$fit_RF, newdata = dat_fold_test, type="prob")[,2], getROC=TRUE)
      }
    }
    list(CVfold)
  }
  
  return (CVfolds)
}





