###
### Brute Force Approach: Test all combinations of proteins
###
doBruteForce <- function(data., predNames) {
  
  allInd <- CreateIndices(length(predNames), length(predNames))
  
  
  require("parallel")
  require("foreach")
  require("doParallel")
  #vignette("parallel")
  
  if (Sys.info()[['sysname']] == "Windows") {
    cl <- makeCluster(detectCores() - 1)
  } else {
    cl <- makeForkCluster(detectCores() - 1)
  }
  registerDoParallel(cl, cores = detectCores() - 1)
  
  allDelta = foreach(i = 1:length(allInd), .packages=c("boot")) %dopar% {
    try({
      fit_Prot <- glm(formula(paste(c("Ectomy_GS >= 7", paste(predNames[allInd[[i]]], collapse="+")), collapse="~")), family=binomial, data=data.)
      cv <- cv.glm(data., fit_Prot)$delta[2]
    })
  }
  stopCluster(cl)
  
  return (allDelta)
}








doFeatureSelectionPlot <- function(RF, GLM, tit="RF Importance of Features", ...) {
  # plot importance of features in RF and selected GLM Features
  colors <- rep("gray", nrow(RF$importance))
  rfnames <- names(rev(sort(RF$importance[,4])))
  colors[which(rfnames %in% names(coef(GLM)))] <- "orange"
  barplot(sort(RF$importance[,4], decreasing = TRUE), main=tit, ylab="Mean Gini Decrease", las=2, cex.names = 0.8, col = colors, ...)
  legend("topright", legend=c("selected by GLM Stepwise Regression"), lty=1, col="orange")
}


doClassifierSurvivalCurves <- function(data, classifier, tit) {
  
  data. <- data
  
  data.$Pred <- predict(classifier, newdata = data.)
  
  if (class(data.$Pred) == "numeric") {
    data.$Pred <- data.$Pred >= 0
  }
  
  Pred.surv <- survfit(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ Pred, data = data.)
  Pred.cox  <-   coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ Pred, data = data.)
  plot(Pred.surv, col = 1:2, ylab="Cumulative Survival", xlab="Months", main=tit, sub=paste(c("n = ", nrow(data.), ",  p=", format(summary(Pred.cox)$logtest[3], digits=3, nsmall=2)), collapse=""))
  if (!is.null(Pred.surv$strata)) {
    legend("bottomleft", legend = names(Pred.surv$strata), col=1:2, lty=1)
  } else if (length(data.$Pred[1])>0) {
    legend("bottomleft", legend = paste(c("Pred=", data.$Pred[1]), collapse = ""), col=1:2, lty=1)
  }
  
  return (summary(Pred.cox)$logtest[3])
}



doPlots <- function(data, predNames, targetname, rseed = 1, doBruteForce = FALSE, testdata = NA) {
  
  data_orig <- data
  data <- removeNas(data, c(predNames,targetname))
  
  ###
  ### RANDOM FOREST
  ###
  require(randomForest)
  
  # Training and OOB-Error on data
  ind_train <- data$Ectomy_GS > 7 | data$Ectomy_GS < 7
  
  set.seed(rseed)
  fit_RF <- randomForest(formula(paste(c("as.factor(", targetname, ")~", paste(predNames, collapse="+")), collapse="")), data.frame(data[ind_train, ]), importance=TRUE, na.action=na.omit, mtry=length(predNames))
  cm <- table(data[ind_train, targetname], fit_RF$votes[,2]>=0.5)
  fit_RF_CV_Error <- sum(diag(cm)) / sum(cm)
  
  ###
  ### LOGISTIC REGRESSION
  ###
  require("boot")
  
  # Training and OOB-Error on data
  fit_GLMFull <- glm(formula(paste(c(targetname, paste(predNames, collapse="+")), collapse="~")), family=binomial, data=data[ind_train, ])
  cv_GLMFull <- my.cv(fit_GLMFull, data=data[ind_train, ])
  cm <- table(data[ind_train, targetname], cv_GLMFull$predictions>=0.5)
  fit_GLMFull_CV_Error <- sum(diag(cm)) / sum(cm)
  
  fit_GLM = step(fit_GLMFull, trace = 0)
  cv_GLM <- my.cv(fit_GLM, data=data[ind_train, ])
  cm <- table(data[ind_train, targetname], cv_GLM$predictions>=0.5)
  fit_GLM_CV_Error <- sum(diag(cm)) / sum(cm)
  
  
  ###
  ### ROC Curve Experiments
  ###
  require(pROC)
  
  
  ###
  ### AUC of PSA alone
  ###
  roc_GLM   <- roc(data[ind_train, targetname], cv_GLM$predictions)
  
  ###
  ### AUC of learned RandomForests
  ###
  roc_RF <- roc(data[ind_train, targetname], fit_RF$votes[,2])
  
  
  ###
  ### Plot the curves
  ###
  
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  old_colors <- palette()
  set.seed(rseed)
  palette(sample(color, 10))
  
  plot(roc_GLM, col=1, main=paste(c("AUC of predictors of ", targetname), collapse=""), sub=paste("n =", sum(ind_train), collapse=""))
  plot(roc_RF, add=TRUE, col=2)
  
  text(0, 0, paste("n =", sum(ind_train), collapse=""),  cex=1.5, pos=2)
  
  
  text(x=0, y=0.52, paste("GLM AUC =", format(roc_GLM$auc, digits=2), collapse=""), pos=2, cex=1.5, col=1)
  text(x=0, y=0.34, paste("RF AUC =", format(roc_RF$auc, digits=2), collapse=""), pos=2, cex=1.5, col=2)
  
  
  
  ###
  ### Survival Plots
  ###
  
  
  #pdf("Pred_Survival2.pdf")
  
  Surv_GLM_train <- doClassifierSurvivalCurves(data, fit_GLM, tit = paste("BR free survival in patients with GLM Predicted", targetname))
  Surv_RF_train <- doClassifierSurvivalCurves(data, fit_RF, tit = paste("BR free survival in patients with RF Predicted", targetname))
  
  doClassifierSurvivalCurves(data[which(data$Ectomy_GS==7), ], fit_GLM, tit = paste("BR free survival in GS7 patients with GLM Predicted", targetname))
  doClassifierSurvivalCurves(data[which(data$Ectomy_GS==7), ], fit_RF, tit = paste("BR free survival in GS7 patients with RF Predicted", targetname))
  
  if (!is.na(testdata) && nrow(testdata)>0) {
    
    Surv_GLM_test <- doClassifierSurvivalCurves(testdata, fit_GLM, tit = paste("BR free survival in TEST patients with GLM Predicted", targetname))
    Surv_RF_test <- doClassifierSurvivalCurves(testdata, fit_RF, tit = paste("BR free survival in TEST patients with RF Predicted", targetname))
    
    doClassifierSurvivalCurves(testdata[which(testdata$Ectomy_GS==7), ], fit_GLM, tit = paste("BR free survival in GS7 TEST patients with GLM Predicted", targetname))
    doClassifierSurvivalCurves(testdata[which(testdata$Ectomy_GS==7), ], fit_RF, tit = paste("BR free survival in GS7 TEST patients with RF Predicted", targetname))
    
    fit_GLMFull_Test <- predict.glm(fit_GLMFull, newdata=testdata, type="response")
    fit_GLM_Test <- predict.glm(fit_GLM, newdata=testdata, type="response")
    fit_RF_Test <- predict(fit_RF, newdata=testdata, type="vote", norm.votes = TRUE)[,2]
    
    cm <- table(testdata[, targetname], fit_GLM_Test>=0.5)
    fit_GLM_Test_Error <- sum(diag(cm)) / sum(cm)
    
    cm <- table(testdata[, targetname], fit_GLMFull_Test>=0.5)
    fit_GLMFull_Test_Error <- sum(diag(cm)) / sum(cm)
    
    cm <- table(testdata[, targetname], fit_RF_Test>=0.5)
    fit_RF_Test_Error <- sum(diag(cm)) / sum(cm)
    
  } else {
    
    Surv_GLM_test <- NA
    Surv_RF_test <- NA
    
    fit_GLMFull_Test <- NA
    fit_RF_Test <- NA
    
    fit_GLM_Test_Error <- NA
    fit_GLMFull_Test_Error <- NA
    fit_RF_Test_Error <- NA
  }
  
  doFeatureSelectionPlot(fit_RF, fit_GLM)
  
  result <- c()
  result$fit_RF <- fit_RF
  result$fit_GLM <- fit_GLM
  result$fit_GLMCV <- cv_GLM
  result$fit_GLMFull <- fit_GLMFull
  result$fit_GLMFullCV <- cv_GLMFull
  result$roc_GLM <- roc_GLM
  result$roc_RF <- roc_RF
  result$Surv_GLM_train <- Surv_GLM_train
  result$Surv_RF_train <- Surv_RF_train
  result$Surv_GLM_test <- Surv_GLM_test
  result$Surv_RF_test <- Surv_RF_test
  result$fit_GLM_test <- fit_GLMFull_Test
  result$fit_RF_Test <- fit_RF_Test
  
  result$fit_GLM_CV_Error <- fit_GLM_CV_Error
  result$fit_GLMFull_CV_Error <- fit_GLMFull_CV_Error
  result$fit_RF_CV_Error <- fit_RF_CV_Error
  result$fit_GLM_Test_Error <- fit_GLM_Test_Error
  result$fit_GLMFull_Test_Error <- fit_GLMFull_Test_Error
  result$fit_RF_Test_Error <- fit_RF_Test_Error
  
  return (result)
}



