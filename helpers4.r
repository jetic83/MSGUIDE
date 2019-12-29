###
### Defines the errormeasure between two vectors y and y_hat (e.g. the accuracy)
###
errormeasure <- function(y, y_hat) {
  # classification error (accuracy)
  cm <- table(y, y_hat)
  e <- sum(diag(cm)) / sum(cm)
  return(e)
}

###
### Defines the AUC errormeasure between two vectors y and y_hat (i.e. the AUC value). y has to be binary, y_hat numerical (probabilities).
###
errormeasureAUC <- function(y, y_hat, doPlot=FALSE, smooth=FALSE, getROC = FALSE, ...) {
  roc_ = roc(y, y_hat, smooth=smooth)
  if (doPlot) {
    plot(roc_, mar=c(5.1,4.1,4.2,2.1), ...)
  }
  if (getROC) {
    return (roc_)
  } else {
    return (auc(roc(y, y_hat)))
  }
}

###
### Defines the Survival errormeasure of a predicted binary outcome y_hat. survivaldata must include
### a column for Time_BCR_free_survival and Status_BCR and be in the same order as y_hat.
###
errormeasureSURV <- function(y_hat, survivaldata) {
  Pred.cox <- tryCatch(coxph(Surv(time=Time_BCR_free_survival, event=Status_BCR, type="right") ~ y_hat, data = survivaldata), error=function(e) {NA})
  if (is.na(Pred.cox)) {
    return (NA)
  } else {
    return (1-summary(Pred.cox)$logtest[3])
  }
}


###
### Removes rows with nas in any predNames column
###
removeNas <- function(data, predNames) {
  nas <- which(is.na(data[,c(predNames)]), arr.ind = TRUE)
  if (length(nas)>0) {
    if (length(dim(nas))==0) {
      data <- data[-unique(nas),]
    } else if(length(nas[,1])>0) {
      data <- data[-unique(nas[,1]),]
    }
  }
  return (data)
}


###
### Finds the timepoint in survival data, after which there are
### as many survivors (event==0) as events (event==1, e.g. death, relapse, etc.)
### If the study ends and there are less events than half of the population,
### the last timepoint is returned.
###
findBalancedTimePoint <- function(eventdata, timedata) {
  dat <- cbind(eventdata, timedata)
  dat <- dat[rowSums(is.na(dat))==0,]
  n <- nrow(dat)
  deaths <- 0
  months <- sort(unique(dat[,2]))
  for (m in months) {
    if (deaths < n/2) {
      deaths <- deaths + sum(dat[which(dat[, 2] == m), 1])
    } else break
  }
  return (m)
}


my.cv <- function(fit, data, classifier="glm") {
  
  predictions <- c()
  
  for (i in 1:nrow(data)) {
    
    if (classifier=="glm") {
      fit <- glm(formula=fit$formula, data=data[-i,], family=fit$family)
      predictions <- c(predictions, predict.glm(fit, newdata=data[i,], type="response"))
    } 
    else if (classifier=="rf") {
      fit <- randomForest(fit$formula, data[-i,], importance=TRUE)
      predictions <- c(predictions, predict(fit, newdata=data[i,], type="vote", norm.votes = TRUE)[,2])
    }
    
  }
  
  list(predictions = predictions)
  
  
}

###
### Helper function for creating all combinations of 1 to p indices. Each index can be 1:n
###
CreateIndices_fast <- function(n, p) {
  # create command
  cmd <- paste(c("sapply(1:(", as.character(n), "-", as.character(p-1), "), function(v", as.character(1), "){PLACEHOLDER})"), collapse="")
  if (p>1) {
    for(i in 2:p) {
      cmd <- sub("PLACEHOLDER", paste(c("sapply((v", as.character(i-1), "+1):(", as.character(n), "-", as.character(p-i), "), function(v", as.character(i), "){PLACEHOLDER})"), collapse=""), cmd)
    }
  }
  cmd <- sub("PLACEHOLDER", paste(c("list(c(", paste(sapply(1:p, function(x){paste(c("v", as.character(x)), collapse="")}), collapse=","), "))"), collapse=""), cmd)
  cmd <- paste(c("unlist(", cmd, ")"), collapse="")
  # evaluate command
  tmp <- eval(parse(text=cmd))
  out <- lapply(1:(length(tmp)/p), function(x){c(tmp[(p*(x-1)+1):(p*x)])})
  return (out)
}


###
### Main function for creating all combinations of 1 to p indices. Each index can be 1:n
###
CreateIndices <- function(n, p) {
  out <- list()
  for (i in 1:p) {
    out <- append(out, CreateIndices_fast(n, i))
  }
  return (out)
}
