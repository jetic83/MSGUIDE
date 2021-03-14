source("helpers.r")


filename <- commandArgs()[5]      
      j  <- as.numeric(commandArgs()[6])
start_i  <- as.numeric(commandArgs()[7])
stop_i   <- as.numeric(commandArgs()[8])
  
# BOOTSTRAPS

GS_RF_PSABX_CVALL <-list() 
GS_RF_PSA_CVALL <-list() 
GS_RF_CVALL <-list() 
PVALS_CVALL <-list() 

for (ind in start_i:stop_i) { 
  
  show(ind)
  predNames <- predNames20[unlist(allInd[[ind]])]
  show(predNames)

  result <- doCV(X, predNames, targetname, rseed = j, k_fold=k, ncores=1)
  
  GS_RF_PSABX_CVALL <- cbind(GS_RF_PSABX_CVALL, unlist(lapply(result, function(x) x$fit_RF_PSABX_error)))
  GS_RF_PSA_CVALL <- cbind(GS_RF_PSA_CVALL, unlist(lapply(result, function(x) x$fit_RF_PSA_error)))
  GS_RF_CVALL <- cbind(GS_RF_CVALL, unlist(lapply(result, function(x) x$fit_RF_error)))
  
  PVALS <- c()
  for (ik in 1:1) {
    PVALS <- c(PVALS, unlist(c(I.pvalue=0, II.pvalue=0, III.pvalue=0, IV.pvalue=0)))
  }
  PVALS_CVALL <- cbind(PVALS_CVALL, unlist(PVALS))
  
}

dimnames(GS_RF_PSABX_CVALL)[[2]] <- as.list(sapply(start_i:stop_i, function(x){paste(predNames20[allInd[[x]]], collapse=", ")}))
write.table(GS_RF_PSABX_CVALL , file=sub("table", "RFCV_PSABX", filename), sep=";", row.names=FALSE)

dimnames(GS_RF_PSA_CVALL)[[2]] <- as.list(sapply(start_i:stop_i, function(x){paste(predNames20[allInd[[x]]], collapse=", ")}))
write.table(GS_RF_PSA_CVALL , file=sub("table", "RFCV_PSA", filename), sep=";", row.names=FALSE)

dimnames(GS_RF_CVALL)[[2]] <- as.list(sapply(start_i:stop_i, function(x){paste(predNames20[allInd[[x]]], collapse=", ")}))
write.table(GS_RF_CVALL , file=sub("table", "RFCV", filename), sep=";", row.names=FALSE)

dimnames(PVALS_CVALL)[[2]] <- as.list(sapply(start_i:stop_i, function(x){paste(predNames20[allInd[[x]]], collapse=", ")}))
write.table(PVALS_CVALL , file=sub("table", "PVALS", filename), sep=";", row.names=FALSE)

