source("helpers.r")

####
####
#### Read the Data
####
alldata = read.csv("alldata_elisa.csv", sep = ",")


predNames20 <- c("PSA_Dens", "BX_GS", "FN1_ELISA", "VTN_ELISA")
predNames20 <- sort(predNames20)


alldata <- removeNas(alldata, unique(c(predNamesRef_PSABX, predNames20, targetname)))

# Prepare CV dataset and Test Dataset

# Train on one batch ProCOC x
train_data_global <- alldata[which(alldata$reading==1 & alldata$ProCOC==batchj), ]
test_data_global <- alldata[which(alldata$reading==1 & alldata$ProCOC!=batchj), ]

require(sampling)
set.seed(rseed)

dat_train <- train_data_global
dat_test  <- test_data_global # for overall validation, will not be used in CV

dat_train <- removeNas(dat_train, unique(c(predNamesRef_PSABX, predNames20, targetname)))
dat_test  <- removeNas(dat_test, unique(c(predNamesRef_PSABX, predNames20, targetname)))

X <- dat_train

# END prepare CV dataset and Test Dataset


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
  
  if (length(unique(dat_fold_test[, targetname]))>1) {
  
    ##
    ## Random Forest
    ##
    require(randomForest)
    set.seed(rseed)
    rf_fit <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(predNamesRef_PSABX, collapse="+")), collapse="")), data.frame(dat_fold_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(predNamesRef_PSABX))
    if (errorfun == "aucerror") {
      RF_refs_PSABX_errors <- c(RF_refs_PSABX_errors, errormeasureAUC(dat_fold_test[, targetname], predict(rf_fit, newdata = dat_fold_test, type="prob")[,2]))
    } 
    RF_refs_PSABX <- c(RF_refs_PSABX, list(rf_fit))
    
    set.seed(rseed)
    rf_fit <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(predNamesRef_PSA, collapse="+")), collapse="")), data.frame(dat_fold_train), importance=TRUE, na.action=na.omit, ntree = ntree, mtry=length(predNamesRef_PSA))
    if (errorfun == "aucerror") {
      RF_refs_PSA_errors <- c(RF_refs_PSA_errors, errormeasureAUC(dat_fold_test[, targetname], predict(rf_fit, newdata = dat_fold_test, type="prob")[,2]))
    } 
    RF_refs_PSA <- c(RF_refs_PSA, list(rf_fit))
  }
}



doClassifierSurvivalFitWithClassifier2 = function (data, classifier, tit, classifierPSA = NA, classifierPSABX = NA, 
                                                   legend = TRUE, doPlot = TRUE, plotSurvivalRatherThanROC = TRUE, 
                                                   plotPSA = TRUE, plotPSABX = TRUE, plotNCCN = TRUE, cutoff = 0.5, 
                                                   cutoffPSA = 0.5, cutoffPSABX = 0.5, ctype = "rf", mark = 3, 
                                                   lwd = 1, ...) 
{
  require(pROC)
  require(survival)
  data. <- data
  if (ctype == "rf") {
    data.$Score <- predict(classifier, newdata = data., type = "prob")[, 
                                                                       2]
  }
  else {
    data.$Score <- predict(classifier, newdata = data., type = "response")
  }
  rocScore = roc(data.[, targetname], data.$Score)
  if (cutoff < 0) {
    cutoff = rocScore$thresholds[which.max(rocScore$specificities + 
                                             rocScore$sensitivities)[1]]
  }
  if (class(data.$Score) == "numeric") {
    data.$Score <- factor(data.$Score >= cutoff)
  }
  levels(data.$Score) <- c("LOW", "HIGH")
  Pred.surv <- survfit(Surv(time = Time_BCR_free_survival, 
                            event = Status_BCR, type = "right") ~ Score, data = data.)
  Pred.cox <- coxph(Surv(time = Time_BCR_free_survival, event = Status_BCR, 
                         type = "right") ~ Score, data = data.)
  if (doPlot) {
    if (plotSurvivalRatherThanROC) {
      plot(Pred.surv, col = 1:2, ylab = "Cumulative Survival", 
           xlab = "Months", main = tit, sub = paste(c("n = ", 
                                                      paste(Pred.surv$n, collapse = " / "), ",  p=", 
                                                      format(summary(Pred.cox)$logtest[3], digits = 3, 
                                                             nsmall = 2)), collapse = ""), lwd = lwd, 
           mark = mark, ...)
      if (legend && !is.null(names(Pred.surv$strata))) {
        leg = names(Pred.surv$strata)
        leglty = c(1, 1)
        legcol = c("black", "red")
        legend("bottomleft", leg, col = legcol, lty = leglty, 
               bty = "n", cex = 0.8)
      }
    }
    else {
      plot(rocScore, col = 1, main = tit, sub = paste(c("n = ", 
                                                        length(rocScore$cases) + length(rocScore$controls), 
                                                        " (p=", format(summary(Pred.cox)$logtest[3], 
                                                                       digits = 3, nsmall = 2), ")"), collapse = ""), 
           mar = c(5.1, 4.1, 4.2, 2.1), ...)
      text(x = 0.4, y = 0.4, labels = paste(c("AUC ", format(auc(rocScore), 
                                                             digits = 3, nsmall = 2)), collapse = ""), pos = 4, 
           col = 1)
    }
    if (plotPSA && !is.na(classifierPSA)) {
      if (ctype == "rf") {
        data.$PSAScore <- predict(classifierPSA, newdata = data., 
                                  type = "prob")[, 2]
      }
      else {
        data.$PSAScore <- predict(classifierPSA, newdata = data., 
                                  type = "response")
      }
      rocPSA = roc(data.[, targetname], data.$PSAScore)
      if (cutoffPSA < 0) {
        cutoffPSA = rocPSA$thresholds[which.max(rocPSA$specificities + 
                                                  rocPSA$sensitivities)[1]]
      }
      if (class(data.$PSAScore) == "numeric") {
        data.$PSAScore <- factor(data.$PSAScore >= cutoffPSA)
      }
      levels(data.$PSAScore) <- c("LOW", "HIGH")
      PredPSA.surv <- survfit(Surv(time = Time_BCR_free_survival, 
                                   event = Status_BCR, type = "right") ~ PSAScore, 
                              data = data.)
      PredPSA.cox <- coxph(Surv(time = Time_BCR_free_survival, 
                                event = Status_BCR, type = "right") ~ PSAScore, 
                           data = data.)
      if (plotSurvivalRatherThanROC) {
        plot(PredPSA.surv, col = 1:2, ylab = "Cumulative Survival", 
             xlab = "Months", main = paste(c("PSA_Dens"), 
                                           collapse = " + "), sub = paste(c("n = ", 
                                                                            paste(PredPSA.surv$n, collapse = " / "), 
                                                                            ",  p=", format(summary(PredPSA.cox)$logtest[3], 
                                                                                            digits = 3, nsmall = 2)), collapse = ""), 
             lwd = lwd, mark = mark, ...)
        if (legend) {
          leg = names(PredPSA.surv$strata)
          leglty = c(1, 1)
          legcol = c("black", "red")
          legend("bottomleft", leg, col = legcol, lty = leglty, 
                 bty = "n", cex = 0.8)
        }
      }
      else {
        plot(rocPSA, col = 1, main = paste(c("PSA_Dens"), 
                                           collapse = " + "), sub = paste(c("n = ", length(rocPSA$cases) + 
                                                                              length(rocPSA$controls), " (p=", format(summary(PredPSA.cox)$logtest[3], 
                                                                                                                      digits = 3, nsmall = 2), ")"), collapse = ""), 
             mar = c(5.1, 4.1, 4.2, 2.1), ...)
        text(x = 0.4, y = 0.3, labels = paste(c("AUC ", 
                                                format(auc(rocPSA), digits = 3, nsmall = 2)), 
                                              collapse = ""), pos = 4, col = "darkgray")
      }
    }
    if (plotPSABX && !is.na(classifierPSABX)) {
      if (ctype == "rf") {
        data.$PSABXScore <- predict(classifierPSABX, 
                                    newdata = data., type = "prob")[, 2]
      }
      else {
        data.$PSABXScore <- predict(classifierPSABX, 
                                    newdata = data., type = "response")
      }
      rocPSABX = roc(data.[, targetname], data.$PSABXScore)
      if (cutoffPSABX < 0) {
        cutoffPSABX = rocPSABX$thresholds[which.max(rocPSABX$specificities + 
                                                      rocPSABX$sensitivities)[1]]
      }
      if (class(data.$PSABXScore) == "numeric") {
        data.$PSABXScore <- factor(data.$PSABXScore >= 
                                     cutoffPSABX)
      }
      levels(data.$PSABXScore) <- c("LOW", "HIGH")
      PredPSABX.surv <- survfit(Surv(time = Time_BCR_free_survival, 
                                     event = Status_BCR, type = "right") ~ PSABXScore, 
                                data = data.)
      PredPSABX.cox <- coxph(Surv(time = Time_BCR_free_survival, 
                                  event = Status_BCR, type = "right") ~ PSABXScore, 
                             data = data.)
      if (plotSurvivalRatherThanROC) {
        plot(PredPSABX.surv, col = 1:2, ylab = "Cumulative Survival", 
             xlab = "Months", main = paste(c("PSA_Dens", 
                                             "BX_GS"), collapse = " + "), sub = paste(c("n = ", 
                                                                                        paste(PredPSABX.surv$n, collapse = " / "), 
                                                                                        ",  p=", format(summary(PredPSABX.cox)$logtest[3], 
                                                                                                        digits = 3, nsmall = 2)), collapse = ""), 
             lwd = lwd, mark = mark, ...)
        if (legend) {
          leg = names(PredPSABX.surv$strata)
          leglty = c(1, 1)
          legcol = c("black", "red")
          legend("bottomleft", leg, col = legcol, lty = leglty, 
                 bty = "n", cex = 0.8)
        }
      }
      else {
        plot(rocPSABX, col = 1, main = paste(c("PSA_Dens", 
                                               "BX_GS"), collapse = " + "), sub = paste(c("n = ", 
                                                                                          length(rocPSABX$cases) + length(rocPSABX$controls), 
                                                                                          " (p=", format(summary(PredPSABX.cox)$logtest[3], 
                                                                                                         digits = 3, nsmall = 2), ")"), collapse = ""), 
             mar = c(5.1, 4.1, 4.2, 2.1), ...)
        text(x = 0.4, y = 0.2, labels = paste(c("AUC ", 
                                                format(auc(rocPSABX), digits = 3, nsmall = 2)), 
                                              collapse = ""), pos = 4, col = "lightblue3")
      }
    }
    if (plotNCCN) {
      rocNCCN = roc(data.[, targetname], data.$NCCN_risk_group_Bx)
      PredNCCN.surv <- survfit(Surv(time = Time_BCR_free_survival, 
                                    event = Status_BCR, type = "right") ~ NCCN_risk_group_Bx > 
                                 2, data = data.)
      PredNCCN.cox <- coxph(Surv(time = Time_BCR_free_survival, 
                                 event = Status_BCR, type = "right") ~ NCCN_risk_group_Bx > 
                              2, data = data.)
      if (plotSurvivalRatherThanROC) {
        plot(PredNCCN.surv, col = 1:2, ylab = "Cumulative Survival", 
             xlab = "Months", main = paste(c("NCCN"), collapse = " + "), 
             sub = paste(c("n = ", paste(PredNCCN.surv$n, 
                                         collapse = " / "), ",  p=", format(summary(PredNCCN.cox)$logtest[3], 
                                                                            digits = 3, nsmall = 2)), collapse = ""), 
             lwd = lwd, mark = mark, ...)
        if (legend) {
          leg = names(PredNCCN.surv$strata)
          leglty = c(1, 1)
          legcol = c("black", "red")
          legend("bottomleft", leg, col = legcol, lty = leglty, 
                 bty = "n", cex = 0.8)
        }
      }
      else {
        plot(rocNCCN, col = 1, main = paste(c("NCCN"), 
                                            collapse = " + "), sub = paste(c("n = ", length(rocNCCN$cases) + 
                                                                               length(rocNCCN$controls), " (p=", format(roc.test(rocScore, 
                                                                                                                                 rocNCCN)$p.value, digits = 3, nsmall = 2), 
                                                                             ")"), collapse = ""), mar = c(5.1, 4.1, 4.2, 
                                                                                                           2.1), ...)
        text(x = 0.4, y = 0.1, labels = paste(c("AUC ", 
                                                format(auc(rocNCCN), digits = 3, nsmall = 2)), 
                                              collapse = ""), pos = 4, col = "palegreen3")
      }
    }
  }
  return(summary(Pred.cox)$logtest[3])
}

doClassifierSurvivalFitWithPrednames2 = function (datatrain, datatest, predNames, targetname, tit = "", 
          ranseed = rseed, ktree = ntree, doPlot = TRUE, cutoff = 0.5, 
          cutoffPSA = 0.5, cutoffPSABX = 0.5, ctype = "rf", plotSurvivalRatherThanROC = TRUE, 
          ...) 
{
  set.seed(ranseed)
  if (ctype == "rf") {
    require(randomForest)
    classifier = randomForest(formula(paste(c("as.factor(", 
                                              targetname, ")~", paste(unique(predNames), collapse = "+")), 
                                            collapse = "")), data.frame(datatrain), importance = TRUE, 
                              na.action = na.omit, ntree = ktree, mtry = length(unique(predNames)))
  }
  set.seed(ranseed)
  if (ctype == "rf") {
    classifierPSA = randomForest(formula(paste(c("as.factor(", 
                                                 targetname, ")~", paste("PSA_Dens", collapse = "+")), 
                                               collapse = "")), data.frame(datatrain), importance = TRUE, 
                                 na.action = na.omit, ntree = ktree)
  }
  set.seed(ranseed)
  if (ctype == "rf") {
    classifierPSABX = randomForest(formula(paste(c("as.factor(", 
                                                   targetname, ")~", paste(unique(c("PSA_Dens", "BX_GS")), 
                                                                           collapse = "+")), collapse = "")), data.frame(datatrain), 
                                   importance = TRUE, na.action = na.omit, ntree = ktree, 
                                   mtry = 2)
  }
  if (is.na(tit) || tit == "") {
    tit = paste(predNames, collapse = " + ")
  }
  return(doClassifierSurvivalFitWithClassifier2(datatest, classifier, 
                                                tit, classifierPSA, classifierPSABX, doPlot = doPlot, 
                                                cutoff = cutoff, cutoffPSA = cutoffPSA, cutoffPSABX = cutoffPSABX, 
                                                ctype = ctype, plotSurvivalRatherThanROC = plotSurvivalRatherThanROC, 
                                                ...))
}



## plots for paper
pdf(paste(c("validate-", paste(predNames20, collapse="-"), ".pdf"), collapse=""))

## Classification error on test_data_global
predNames    <- predNames20

# classifier retrained on whole training set


# best classifier from cross validation
result = doCV(X, predNames, targetname, rseed = rseed, k_fold=k, ncores=1)
GS_RF_CVAll_ = unlist(lapply(result, function(x) x$fit_RF_error))
inds_fold <- train_inds[[which.max(GS_RF_CVAll_)]]
inds_oob_fold <- test_inds[[which.max(GS_RF_CVAll_)]]
dat_fold_train <- X[inds_fold, ]
dat_fold_test  <- X[inds_oob_fold, ]
set.seed(rseed)
if (ctype=="rf") {
  val_RF = randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(unique(c(predNames)), collapse="+")), collapse="")), X[inds_fold, ], importance=TRUE, na.action=na.omit, mtry=length(unique(c(predNames))))
} 

if (errorfun == "aucerror") {
  if (ctype=="rf") {
    val_RF_error <- errormeasureAUC(test_data_global[, targetname], predict(val_RF, newdata = test_data_global, type="prob")[,2])
  }
} 

cutoffs <- comparisonAUCPlotsWReferences(unique(c(predNames)), X[inds_fold, ], test_data_global, r_seed = rseed, ctype=ctype)
cutoff <- cutoffs[1]
cutoffPSA <- cutoffs[2]
cutoffPSABX <- cutoffs[3]


# do the AUC on the 118 subset patients only
train_data_global$usable = !is.na(train_data_global$Time_BCR_free_survival) & !is.na(train_data_global$PSA_Dens)
doClassifierSurvivalFitWithPrednames(X[inds_fold, ], train_data_global[train_data_global$usable,], predNames, targetname, paste(c(paste(predNames, collapse=" + ")), collapse=" "), ctype=ctype, plotSurvivalRatherThanROC=FALSE)

## Survival Curve on test_data_global
doClassifierSurvivalFitWithPrednames(X[inds_fold, ], test_data_global, predNames, targetname, paste(c(paste(predNames, collapse=" + ")), collapse=" "), cutoff=cutoff, cutoffPSA=cutoffPSA, cutoffPSABX=cutoffPSABX, ctype=ctype, lwd=4)
doClassifierSurvivalFitWithPrednames2(X[inds_fold, ], test_data_global, predNames, targetname, paste(c(paste(predNames, collapse=" + ")), collapse=" "), cutoff=cutoff, cutoffPSA=cutoffPSA, cutoffPSABX=cutoffPSABX, ctype=ctype, lwd=4)
doClassifierSurvivalFitWithPrednames(X[inds_fold, ], test_data_global, predNames, targetname, paste(c(paste(predNames, collapse=" + ")), collapse=" "), cutoff=cutoff, cutoffPSA=cutoffPSA, cutoffPSABX=cutoffPSABX, ctype=ctype, plotSurvivalRatherThanROC=FALSE)

dev.off()



