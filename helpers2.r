require(randomForest)
require(pROC)
require(MASS)
require(igraph)
require(irr)
require(rms)
require(AER)
require(entropy)




################################################################################
##
##  Helper Functions
##
##  Peter J. Sch_ffler
##  ETH Zurich, Institute for Computational Science
##  peter.schueffler@inf.ethz.ch
##
################################################################################


################################################################
### Brutus Jobs processing
################################################################
# Checks if no jobs are pending anymore and more than skip jobs are running.
jobsnmax <- function() {
  defnmax <- 288
  cmd <- "busers"
  out <- system(cmd, intern = TRUE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
  if (grepl("Cannot connect to LSF. Please wait ...", out[1])) {
    return(defnmax)
  }
  strsplit(out[1], " +")[[1]]
  out <- strsplit(out[2], " +")[[1]]
  return (tryCatch(as.numeric(out[3]), warning=function(e){defnmax}, error=function(e){defnmax}))
}

jobsFinished <- function(skip=0) {
  cmd <- "busers"
  out <- system(cmd, intern = TRUE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
  if (grepl("Cannot connect to LSF. Please wait ...", out[1])) {
    return(FALSE)
  }
  strsplit(out[1], " +")[[1]]
  out <- strsplit(out[2], " +")[[1]]
  return (out[5]=="0" && !is.na(as.numeric(out[6])) && as.numeric(out[6])<=skip)
}

waitUntilJobsFinished <- function(skip=0) {
  Sys.sleep(5)
  while(!jobsFinished(skip)) {
    Sys.sleep(60)
  }
  return (TRUE);
}

informUserViaEmail <- function(address = "peter.schueffler@inf.ethz.ch") {
  #require(source("http://rtm.wustl.edu/code/sendEmail.R")
  # make sure that field subject is only one string without any spaces!!! (Or use \"\" as substring with spaces)
  #sendEmail(subject = "\"Euler Job Finished\"", text = paste("Your Euler jobs have finished on", Sys.time(), "."), address = address)
  require(sendmailR)
  sendmail(from="<peschuef@brutus.ethz.ch>", to=address, subject = "Euler Job Finished", msg = paste("Your Euler job has finished on", Sys.time(), "."))
}


################################################################
### Plot error bars
################################################################
error.bar <- function(x, y, upper, lower=upper, length=0.1, horiz=FALSE, ...) {
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  if (horiz) {
    arrows(y+upper,x, y-lower, x, angle=90, code=3, length=length, ...)
  } else {
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }
}


################################################################
### Plot advanced boxplots
################################################################
myboxplot <- function(x, ci=0.95, ...) {
  out <- boxplot(x, ...)
  if (typeof(x)!="list") {
    x <- list(x)
  }
  means <- sapply(x, function(y) mean(y, na.rm=TRUE))
  errors <- sapply(x, function(y) qt(ci,df=length(y[which(!is.na(y))])-1)*sd(y, na.rm=TRUE)/sqrt(length(y[which(!is.na(y))])))
  segments(1:length(x),means-errors,1:length(x),means+errors, lwd=10, col="lightgrey")
  points(means,col="black", cex=2, pch="+")
  #text(1:length(x), min(unlist(x), na.rm=TRUE), labels=out$n, cex=1.5)
  return (out)
}


################################################################
### Cross Validation
################################################################

loocv <- function(samples, target, mformula, regression = TRUE)
{
  samplesize <- dim(samples)[1]
  errors <- c()
  targetindex <- which(names(samples)==target)
  
  for (i in 1:samplesize)
  {
    # Split the training- / testsets
    trains <- samples[-i,]
    tests <- samples[i,]
    
    # Train the model on trainset
    submodel <- glm(mformula, data=trains)
    
    # Classify the test sets and store the error
    # if regression: the error is the RMS error
    # if classification: the error is the amount of misclassified testsamples
    error <- 0
    for (j in 1:dim(tests)[1])
    {
      testerror <- 0
      if (regression)
      {
        testerror <- (predict.glm(submodel, tests[j,]) - tests[j,targetindex])^2
      } else {
        testerror <- (predict.glm(submodel, tests[j,]) - tests[j,targetindex])^2
      }
      
      error <- error + testerror
    }
    if (regression)
    {
      errors[i] <- sqrt(error/dim(tests)[1])
    } else {
      errors[i] <- sqrt(error/dim(tests)[1])
    }
  }
  
  return(errors)
}


# Gives the first n proteins out of this list. Since proteins are measured several times,
# it's not only proteins[1:20]. Proteins occur only 1 time in the return value.
# If exactName is false, only the protein-name is considered (i.e. "ICAM1_MRM3_NP" and "ICAM1_MRM3_NPr"
# will lead to "ICAM1"). If it's true, the example will lead to "ICAM1_MRM3_NP". This matters if you
# want to compare two firstNProteins-lists, depending on which of "ICAM1_MRM3_NP" and "ICAM1_MRM3_NPr"
# occurs first. proteins is a vector of size >=n of the protein names.
firstNProteins <- function(proteins, n=20, exactName = FALSE) {
  ind=which(proteins=="(Other)")
  if (length(ind)>0) {proteins=proteins[-ind]}
  proteinnames <- sapply(proteins, function(x){strsplit(x,"_")[[1]][1]})                                 
  singles <- which(!duplicated(proteinnames))
  if (n > length(singles)) stop("Error in firstNProteins: length(proteins) < n")
  if (exactName) proteins[singles[1:n]]
  else proteinnames[singles[1:n]]
}

# Returns the first exactname of exactnamelist that corresponds to proteinname.
proteinname2Exactname <- function(proteinname, exactnamelist) {
  if (proteinname %in% exactnamelist) proteinname
  else {
    proteinnames <- sapply(exactnamelist, function(x){strsplit(x,"_")[[1]][1]})
    exactnamelist[which(proteinnames == proteinname)[1]]
  }
}

cname <- function(proteinname) {
  proteinname <- gsub("_MRM2", "", proteinname) 
  proteinname <- gsub("_MRM3", "", proteinname)
  proteinname
}

removeBadColumns <- function(X, t_correlation = 0.9, cormethod="pearson") {
  remcol <- c(dim(X)[2]+1)
  for (i in 1:(dim(X)[2]-1)) {
    if (!(i %in% remcol)) {
      for (j in (i+1):dim(X)[2]) {
        v <- abs(cor(X[,i], X[,j], use="pairwise.complete.obs", method=cormethod))
        if (is.na(v)) { v <- 0 }
        if (v > t_correlation) {
          remcol <- c(remcol, j)
        }
      }
    }
  }
  X_corrected <- X[, -unique(remcol)]
  
  X_corrected <- X_corrected[,which(!grepl("_s_", names(X_corrected)))]
  X_corrected <- X_corrected[,which(!grepl("_FT_", names(X_corrected)))]
  X_corrected <- X_corrected[,which(!grepl("_BT_", names(X_corrected)))]
  return(X_corrected)
}

markers2index <- function(predNames) {
  predNames <- sort(predNames)
  ind <- 0
  ind <- which(unlist(lapply(bsindAll, function(x){length(predNames20[x])==length(predNames) && sum(!(predNames20[x] %in% predNames))==0})))
  ind
}


# Einige Marker, mit Bootstrap validiert
valBS <- function(predNames = c("CRP", "Leukocytosis"), X) {
  indRowNA <- which(apply(is.na(d.f(X[,predNames])),1,sum) > 0)
  if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
  k = 100 # Bootsraps k times
  n = dim(X[-indRowNA,])[1]
  allind <- c(1:n)
  errors <- c()
  for (i in 1:k) {
    # generiere Trainingsset
    bsind <- sample(allind, n, replace=TRUE)
    
    #generiere Testset
    oob <- allind[is.na(match(allind, bsind))]
    testsamples <- X[-indRowNA,][oob,]
    
    T <- as.factor(X[-indRowNA,]$Dx[bsind])
    Y <- d.f(d.f(X[-indRowNA,predNames])[bsind,])
    names(Y) <- predNames
    bsfit <- glm(T ~ ., data=d.f(Y), family=binomial)
    # klassifiziere und speichere den error
    result <- factor(predict(bsfit, newdata=testsamples, type="response") >= 0.5)
    levels(result)[1]<-"BPH"
    levels(result)[2]<-"Cancer"
    error <- 1- sum(result == testsamples$Dx) / dim(testsamples)[1]
    # alle 100 BS-Errors f_r alle 100 Modelle werden gespeichert
    errors <- append(errors, error)
  }
  boxplot(errors, main=paste(c("Bootstrap Error of Linear Model for Dx\n", predNames), collapse=" "), ylab="Bootstrap Error (0-1)", xlab=paste(c("Median =", median(errors)), collapse=" "))
}



# Adds two performances (x-values are merged, y values are added accordingly)
get.perf.value.at <- function(x, perf) {
  x1 <- slot(perf, "x.values")[[1]]
  y1 <- slot(perf, "y.values")[[1]]
  
  m <- 0
  if (is.element(x, x1)) {
    m <- max(y1[which(x1==x)])
  } else {
    ya <- x1[rev(which(x1<x))[1]]
    yb <- x1[(which(x1>x))[1]]
    a <- y1[rev(which(x1<x))[1]]
    b <- y1[(which(x1>x))[1]] 
    
    m <- a + (b-a)*(x-ya)/(yb-ya)
  }
  return(m)
}  

# average performances (x-values are merged, y values are averaged accordingly)
avg.performances <- function(perfs) {
  perf <- perfs[[1]]
  
  xs <- c()
  for (i in 1:length(perfs)) {
    xs <- c(xs, slot(perfs[[i]], "x.values")[[1]])
  }
  
  xs <- sort(unique(xs))
  ys <- xs-xs
  for (i in 1:length(xs)) {
    ys[i] <- median(sapply(1:length(perfs), function(x) {get.perf.value.at(xs[i], perfs[[x]])}))
  }
  
  # form the curve xs/ys to a step-curve
  xs <- c(unlist(lapply(xs, function(x) {c(x, x)})))
  ys <- c(unlist(lapply(ys, function(x) {c(x, x)})))
  ys <- c(0, ys[1:length(ys)-1])
  
  slot(perf, "x.values")[[1]] <- xs
  slot(perf, "y.values")[[1]] <- ys
  
  return(perf)
}

# Adds two rocs (x-values are merged, y values are added accordingly)
get.sens.at <- function(x, ro) {
  x1 <- ro$specificities
  y1 <- ro$sensitivities
  
  m <- 0
  if (is.element(x, x1)) {
    m <- max(y1[which(x1==x)])
  } else {
    ya <- x1[rev(which(x1<x))[1]]
    yb <- x1[(which(x1>x))[1]]
    a <- y1[rev(which(x1<x))[1]]
    b <- y1[(which(x1>x))[1]] 
    
    m <- a + (b-a)*(x-ya)/(yb-ya)
  }
  return(m)
}  

# average roc (x-values are merged, y values are averaged accordingly)
avg.roc <- function(rocs) {
  
  # filter NAs (non-rocs)
  rocs <- rocs[which(!is.na(rocs))]
  
  ro <- rocs[[1]]
  
  #xs <- c()
  #for (i in 1:length(rocs)) {
  #  xs <- c(xs, rocs[[i]]$specificities)
  #}
  #xs <- sort(unique(xs))
  
  xs <- (0:100)/100
  ys <- xs-xs
  for (i in 1:length(xs)) {
    ys[i] <- median(sapply(1:length(rocs), function(x) {get.sens.at(xs[i], rocs[[x]])}))
  }
  
  ys[1] = 0
  ys[length(xs)] = 1
  
  ## form the curve xs/ys to a step-curve
  #xs <- c(unlist(lapply(xs, function(x) {c(x, x)})))
  #ys <- c(unlist(lapply(ys, function(x) {c(x, x)})))
  #ys <- c(0, ys[1:length(ys)-1])
  
  ro$specificities <- xs
  ro$sensitivities <- ys
  
  
  # recalculate auc
  ro$auc = 0
  for (i in 2:length(xs)) {
    ro$auc = ro$auc + (1/length(xs) * (ys[i]+ ys[i-1]) / 2)
  }
  
  return(ro)
}


# Histogramm ueber Praediktoren
predHist <- function(oneToX = 1000, tit = "", orderLike = oneToX, short.names = FALSE, cex.names=0.6, mai=c( 1.02, 1.35, 0.82, 0.42), conditioned=c(), conditioned_inclusive=TRUE, fromX=1, doPlot=TRUE, normalized=FALSE, barcol="orange", ...) {  
  if (!is.na(orderLike)) {
    orderpreds <- unlist(lapply(bsindAll[fromX:orderLike], function(x){predNames20[unlist(x)]}))
    his <- predNames20
    his <- sapply(his, function(x){sum((orderpreds==x))})
    orderpreds <- (names(sort(his)))
  }
  
  allpreds <- lapply(bsindAll[fromX:oneToX], function(x){
    if (length(conditioned)==0) {
      return (predNames20[unlist(x)])
    } else {
      if (conditioned_inclusive) {
        if (all(conditioned %in% predNames20[unlist(x)])) {
          return (predNames20[unlist(x)])
        } else {
          return (c())
        }
      } else {
        if (any(conditioned %in% predNames20[unlist(x)])) {
          return (c())
        } else {
          return (predNames20[unlist(x)])
        }
      }
    }
  })                                                                                                                                                                                         
  if (length(conditioned)>0) {
    n <- sum(sapply(allpreds, function(x){as.numeric(length(x)>0)}))
  } else {
    n <- oneToX-fromX+1
  }
  allpreds <- unlist(allpreds)
  
  his <- predNames20
  his <- sapply(his, function(x){sum((allpreds==x))/n})
  if (normalized) {
    his <- his/sum(his)
  }
  if (!is.na(orderLike)) his <- his[orderpreds]
  
  if (short.names) {
    names(his) <- shortenNames(names(his), reallyshort=TRUE)
  }
  
  if (doPlot) {
    par(mai=mai)
    subt <- ""
    if (length(conditioned)>0) {
      cond <- "include"
      if (!conditioned_inclusive) cond <- "exclude" 
      subt <- paste(c("in models which ", cond, " ", paste(c(conditioned), collapse=", "), ", n=", as.character(n)), collapse="")
    }
    barplot(100*his, horiz=TRUE, col=barcol, las=1, cex.names=cex.names, main=tit, sub=subt, ...)
    par(mai=c(1.02, 0.82, 0.82, 0.42))
  }
  return(his)
}

# Histogramm ueber Praediktoren
predHist2 <- function(allpreds, tit = "", short.names = FALSE, cex.names=0.6, mai=c( 1.02, 1.35, 0.82, 0.42), n=-1, normalized=FALSE, doPlot=TRUE, barcol="orange", ...) {
  his <- unique(unlist(allpreds))
  his <- sapply(his, function(x){sum(unlist(allpreds)==x)})
  his <- sort(his)
  
  if(n>0) {
    his <- his/n
  }
  if (normalized) {
    his <- his/sum(his)
  }
  
  if (short.names) {
    names(his) <- shortenNames(names(his), reallyshort=TRUE)
  }
  
  if (doPlot) {
    par(mai=mai)
    if (n>0) {
      his <- 100*his
    }
    barplot(his, horiz=TRUE, col=barcol, las=1, cex.names=cex.names, main=tit, ...)
    par(mai=c(1.02, 0.82, 0.82, 0.42))
  }
  return(his)
}

# Histogramm ueber Praediktorenpaare
predHist3 <- function(oneToX = 1000, tit = "", orderLike = oneToX, cex.names=0.6, mai=c( 1.02, 1.35, 0.82, 0.42), fromX=1, doPlot=TRUE, normalized=FALSE, plotOneToX=-1, barcol="orange", ...) {
  if (!is.na(orderLike)) {
    orderpreds <- unlist(lapply(bsindAll[fromX:orderLike], function(x){return (if (length(x)<2) NA else apply(combn(predNames20[x], 2), 2, function(x){paste(x, collapse=",")}))}))
    if (length(which(is.na(orderpreds)))>0) orderpreds <- orderpreds[-which(is.na(orderpreds))]
    his <- apply(combn(length(predNames20), 2), 2, function(x){paste(predNames20[x], collapse=",")})
    his <- sapply(his, function(x){sum((orderpreds==x))})
    orderpreds <- (names(sort(his)))
  }
  
  allpreds <- unlist(lapply(bsindAll[fromX:oneToX], function(x){return (if (length(x)<2) NA else apply(combn(predNames20[x], 2), 2, function(x){paste(x, collapse=",")}))}))
  if (length(which(is.na(allpreds)))>0) allpreds <- allpreds[-which(is.na(allpreds))]
  n <- oneToX-fromX+1
  
  his <- apply(combn(length(predNames20), 2), 2, function(x){paste(predNames20[x], collapse=",")})
  his <- sapply(his, function(x){sum((allpreds==x))/n})
  if (normalized) {
    his <- his/sum(his)
  }
  if (!is.na(orderLike)) his <- his[orderpreds]
  
  if (doPlot) {
    par(mai=mai)
    subt <- ""
    if (plotOneToX<1) {
      plotOneToX <- length(his)
    }
    barplot(100*his[(length(his)-plotOneToX+1):length(his)], horiz=TRUE, col=barcol, las=1, cex.names=cex.names, main=tit, sub=subt, ...)
    par(mai=c(1.02, 0.82, 0.82, 0.42))
  }
  return(his)
}

# Graph for independence of features
predHist4 <- function(oneToX = 1000, tit = "", cex.names=0.6, mai=c( 1.02, 1.35, 0.82, 0.42), fromX=1, th=0.00, plot_independence=FALSE, normalized=FALSE, doPlot=TRUE, ...) {
  ph1 <- predHist(oneToX=oneToX, tit=tit, orderLike=NA, cex.names=cex.names, mai=mai, fromX=fromX, doPlot=FALSE, normalized=normalized, ...)
  ph2 <- predHist3(oneToX=oneToX, tit=tit, orderLike=NA, cex.names=cex.names, mai=mai, fromX=fromX, doPlot=FALSE, normalized=normalized, ...)
  
  vertices <- data.frame(name = names(ph1), P = ph1, stringsAsFactors = TRUE)
  
  edges.string <- sapply(names(ph2),function(x){strsplit(x, ",")}, USE.NAMES=FALSE)
  edges <- data.frame(from=sapply(edges.string, function(x){x[1]}),
                      to = sapply(edges.string, function(x){x[2]}),
                      Pij = ph2, 
                      stringsAsFactors = TRUE)
  
  if (plot_independence) {
    ## for Pij < Pi*Pj
    edges$Weight <- apply(edges, 1, function(x){as.numeric(x[3]) - (ph1[x[1]] * ph1[x[2]])})   
  } else {
    ## for log scale
    edges$Weight <- apply(edges, 1, function(x){ifelse(as.numeric(x[3])==0, 0, (as.numeric(x[3]))*log2((as.numeric(x[3]))/(ph1[x[1]] * ph1[x[2]])))})
    vertices$P <- sapply(vertices$P, function(x){-(x)*log2((x))})   
  }
  edges$Weight[which(is.na(edges$Weight))] <- 0
  vertices$P[which(is.na(vertices$P))] <- 0
  
  if (th>0) {
    ind <- which(abs(edges$Weight)<th)
    if (length(ind>0)) {
      edges <- edges[-ind,]
    }
  }
  ind <- which(edges$Weight==0)
  if (length(ind>0)) {
    edges <- edges[-ind,]
  }
  
  row.names(edges) <- NULL
  row.names(edges) <- NULL
  g <- graph.data.frame(edges, directed=FALSE, vertices=vertices)
  E(g)$color <- 2+sign(edges$Weight)
  if (plot_independence) {
    E(g)$width <- 400*abs(edges$Weight)
    V(g)$size <- 50*vertices$P
  } else {
    if (normalized) {
      E(g)$width <- 300*abs(edges$Weight)
      V(g)$size <- 100*vertices$P
    } else {
      E(g)$width <- 500*abs(edges$Weight)
      V(g)$size <- 100*vertices$P
    }
  }
  if (doPlot) {
    op <- par(no.readonly = TRUE)
    set.seed(1)
    plot(g, layout=layout.circle)
    title(tit, sub=paste(c("Threshold z = ", th), collapse=""))
    if (plot_independence) {
      ## For Pij < Pi*pj
      if (normalized) {
        legend("bottomleft", inset=c(-0.03,-0.1), legend=c(expression(paste("p(F"[i], ",F"[j], ") - p(F"[i], ") * p(F"[j], ")", sep="")), expression(paste("Vertice Size = p(F"[i], ")", sep=""))), lty=1, col=c(3,0))
      } else {
        legend("bottomleft", inset=c(-0.03,-0.1), legend=c(expression(paste("p(F"[i], ",F"[j], ") < p(F"[i], ") * p(F"[j], ")", sep="")), expression(paste("p(F"[i], ",F"[j], ") > p(F"[i], ") * p(F"[j], ")", sep="")), expression(paste("Vertice Size = p(F"[i], ")", sep=""))), lty=1, col=c(1,3,0))
      }
    } else {
      ## For log scale
      if (normalized) {
        legend("bottomleft", inset=c(-0.03,-0.1), legend=c(expression(paste("p(F"[i], ",F"[j], ") * log(p(F"[i], ",F"[j], ") / (p(F"[i], ") * p(F"[j], ")))", sep="")), expression(paste("Vertice Size = -p(F"[i], ") * log(p(F"[i], "))", sep=""))), lty=1, col=c(3,0))
      } else {
        legend("bottomleft", inset=c(-0.03,-0.13), legend=c(expression(paste("p(F"[i], ",F"[j], ") * log(p(F"[i], ",F"[j], ") / (p(F"[i], ") * p(F"[j], "))) < 0", sep="")), expression(paste("p(F"[i], ",F"[j], ") * log(p(F"[i], ",F"[j], ") / (p(F"[i], ") * p(F"[j], "))) > 0", sep="")), expression(paste("Vertice Size = -p(F"[i], ") * log(p(F"[i], "))", sep=""))), lty=1, col=c(1,3,0))
      }
    }
    par(op)
  }
  return(g)
}

predHist4Matrix <- function(oneToX = 1000, tit = "", cex.names=0.6, mai=c( 1.02, 1.35, 0.82, 0.42), fromX=1, th=0.00, plot_independence=FALSE, normalized=FALSE, doPlot=TRUE, ...) {
  ph1 <- predHist(oneToX=oneToX, tit=tit, orderLike=NA, cex.names=cex.names, mai=mai, fromX=fromX, doPlot=FALSE, normalized=normalized, ...)
  ph2 <- predHist3(oneToX=oneToX, tit=tit, orderLike=NA, cex.names=cex.names, mai=mai, fromX=fromX, doPlot=FALSE, normalized=normalized, ...)
  
  M <- matrix(0.0, length(ph1), length(ph1), dimnames=list(names(ph1), names(ph1)))
  
  vertices <- data.frame(name = names(ph1), P = ph1, stringsAsFactors = TRUE)
  
  edges.string <- sapply(names(ph2),function(x){strsplit(x, ",")}, USE.NAMES=FALSE)
  edges <- data.frame(from=sapply(edges.string, function(x){x[1]}),
                      to = sapply(edges.string, function(x){x[2]}),
                      Pij = ph2, 
                      stringsAsFactors = TRUE)
  
  if (plot_independence) {
    ## for Pij < Pi*Pj
    edges$Weight <- apply(edges, 1, function(x){as.numeric(x[3]) - (ph1[x[1]] * ph1[x[2]])})   
  } else {
    ## for log scale
    edges$Weight <- apply(edges, 1, function(x){ifelse(as.numeric(x[3])==0, 0, (as.numeric(x[3]))*log2((as.numeric(x[3]))/(ph1[x[1]] * ph1[x[2]])))})
    vertices$P <- sapply(vertices$P, function(x){-(x)*log2((x))})   
  }
  edges$Weight[which(is.na(edges$Weight))] <- 0
  vertices$P[which(is.na(vertices$P))] <- 0
  
  if (nrow(edges)>0) {
    ind <- which(edges$Weight==0)
    if (length(ind>0)) {
      edges <- edges[-ind,]
    }
    
    for (i in 1:nrow(edges)) {
      M[as.character(edges[i,2]), as.character(edges[i,1])] <- edges[i,4]
    }
  }
  M[upper.tri(M)] <- NA
  diag(M) <- NA
  
  if (doPlot) {
    op <- par(no.readonly = TRUE)
    heatmap(M, Rowv=NA, Colv=NA, scale="none", cexRow=1.4, cexCol=1.4, margins=c(15,15))
    image.plot(legend.only=TRUE, zlim= c(min(M, na.rm=TRUE),max(M, na.rm=TRUE)), col=heat.colors(12), smallplot=c(0.5,0.54,0.35,0.71))
    
    #image.plot(legend.only=TRUE, zlim= c(min(M, na.rm=TRUE),max(M, na.rm=TRUE)), col=heat.colors(12), smallplot=c(0.8,0.84,0.02,0.3))
    #text(x=0.4, y=0.3, labels=expression(paste("p(F"[i], ",F"[j], ") - p(F"[i], ") * p(F"[j], ") = ", sep="")), pos=4, srt=90, cex=1.4)
    
    if (FALSE) {
      if (plot_independence) {
        ## For Pij < Pi*pj
        if (normalized) {
          legend("bottomright", legend=c(expression(paste("p(F"[i], ",F"[j], ") - p(F"[i], ") * p(F"[j], ")", sep=""))), fill=heat.colors(12)[12])
        } else {
          legend("bottomright", legend=c(expression(paste("p(F"[i], ",F"[j], ") > p(F"[i], ") * p(F"[j], ")", sep="")), expression(paste("p(F"[i], ",F"[j], ") < p(F"[i], ") * p(F"[j], ")", sep=""))), fill=heat.colors(12)[c(12,1)])
        }
      } else {
        ## For log scale
        if (normalized) {
          legend("bottomright", legend=c(expression(paste("p(F"[i], ",F"[j], ") * log(p(F"[i], ",F"[j], ") / (p(F"[i], ") * p(F"[j], ")))", sep=""))), fill=heat.colors(12)[12])
        } else {
          legend("bottomright", legend=c(expression(paste("p(F"[i], ",F"[j], ") * log(p(F"[i], ",F"[j], ") / (p(F"[i], ") * p(F"[j], "))) > 0", sep="")), expression(paste("p(F"[i], ",F"[j], ") * log(p(F"[i], ",F"[j], ") / (p(F"[i], ") * p(F"[j], "))) < 0", sep=""))), fill=heat.colors(12)[c(12,1)])
        }
      }
    }
    par(op)
  }
  return(M)
}

# Confusion Table of independence of Features
# P(Fi=1, Fj=1)* P(Fi=1, Fj=0)
# P(Fi=0, Fj=1)   P(Fi=0, Fj=0)
# * = relative number of classifiers with feature i and with feature j.
# Iij = SUM(Fi in {0,1}) SUM(Fj in {0,1}) [ P(Fi, Fj) * log( P(Fi, Fj) / P(Fi)*P(Fj) ) ]
# Draw a Graph with Iij as edge width and nodes radius -P(i)*log(Pi)
predHist5 <- function(oneToX = 1000, tit = "", cex.names=0.6, mai=c( 1.02, 1.35, 0.82, 0.42), fromX=1, th=0.00, doPlot=TRUE, ...) {
  Fi <- predHist(oneToX=oneToX, tit=tit, orderLike=NA, cex.names=cex.names, mai=mai, fromX=fromX, doPlot=FALSE, normalized=FALSE, ...)
  Fi1Fj1 <- predHist3(oneToX=oneToX, tit=tit, orderLike=NA, cex.names=cex.names, mai=mai, fromX=fromX, doPlot=FALSE, normalized=FALSE, ...)
  
  #Fi is P(Fi=1)
  #Fi1Fj1 is P(Fi=1, Fj=1)
  # calculate Fi1Fj0, Fi0Fj1, Fi0Fj0 and Iij
  Fi1Fj0 <- Fi1Fj1
  Fi0Fj1 <- Fi1Fj1
  Fi0Fj0 <- Fi1Fj1
  Iij <- Fi1Fj1
  for (i in 1:length(Fi1Fj1)) {
    pFi <- Fi[strsplit(names(Fi1Fj1[i]), ",")[[1]][1]]
    pFj <- Fi[strsplit(names(Fi1Fj1[i]), ",")[[1]][2]]
    Fi1Fj0[i] <- pFi - Fi1Fj1[i]
    Fi0Fj1[i] <- pFj - Fi1Fj1[i]
    Fi0Fj0[i] <- 1 - Fi1Fj1[i] - Fi1Fj0[i] - Fi0Fj1[i]
    cm <- matrix(c(Fi1Fj1[i], Fi1Fj0[i], Fi0Fj1[i], Fi0Fj0[i]), ncol=2)
    Iij[i] <- mutual_information_from_table(cm, unit="log2")$mi_normalized
  }
  Iij[which(is.na(Iij))] <- 0
  
  vertices <- data.frame(name = names(Fi), P = Fi, stringsAsFactors = TRUE)
  vertices$P <- sapply(vertices$P, function(x){-x*log(x)})
  vertices$P[which(is.na(vertices$P))] <- 0
  
  edges.string <- sapply(names(Fi1Fj1),function(x){strsplit(x, ",")}, USE.NAMES=FALSE)
  edges <- data.frame(from=sapply(edges.string, function(x){x[1]}),
                      to = sapply(edges.string, function(x){x[2]}),
                      Pij = Iij, 
                      stringsAsFactors = TRUE)
  edges$Weight <- Iij
  
  if (th>0) {
    ind <- which(abs(edges$Weight)<th)
    if (length(ind>0)) {
      edges <- edges[-ind,]
    }
  }
  ind <- which(edges$Weight==0)
  if (length(ind>0)) {
    edges <- edges[-ind,]
  }
  
  if (doPlot) {
    op <- par(no.readonly = TRUE)
    g <- graph.data.frame(edges, directed=FALSE, vertices=vertices)
    E(g)$color <- 2+sign(edges$Weight)
    E(g)$width <- 300*abs(edges$Weight)
    V(g)$size <- 100*vertices$P
    set.seed(1)
    plot(g, layout=layout.circle)
    title(tit, sub=paste(c("Threshold z = ", th), collapse=""))
    legend("bottomleft", inset=c(0,-0.1), legend=c(expression(paste("I"["ij"], " = ", Sigma[paste("F"[i], ",F"[j], " in {0,1}", sep="")], "p(F"[i], ",F"[j], ") * log(p(F"[i], ",F"[j], ") / p(F"[i], ")*p(F"[j], "))", sep="")), expression(paste("Vertice Size = -p(F"[i], ") * log(p(F"[i], "))", sep=""))), lty=1, col=c(3,0))
    par(op)
  }
  
  return(g)
  
}

preHist6 <- function (predNames, modelsAll, whichModels = NA, cex.names = 0.6, 
                      mai = c(1.02, 1.35, 0.82, 0.42), doPlot = TRUE, tit = "", 
                      barcol = "orange", ...) 
{
  if (all(is.na(whichModels))) {
    whichModels = 1:length(modelsAll)
  }
  predhists = c()
  for (i in 1:length(predNames)) {
    predhists = c(predhists, list(c()))
  }
  names(predhists) = predNames
  for (modInd in whichModels) {
    model = modelsAll[[modInd]]
    correlation = median(model$cv, na.rm = TRUE)
    for (pred in model$predNames) {
      predhists[[pred]] = c(predhists[[pred]], correlation)
    }
  }
  predhists = predhists[order(sapply(predhists, function(x) {
    mean(x, na.rm = TRUE)
  }), na.last = FALSE)]
  predhist = predhists
  predhist$mean = unlist(sapply(predhists, function(x) {
    mean(x, na.rm = TRUE)
  }))
  predhist$sd = unlist(sapply(predhists, function(x) {
    sd(x, na.rm = TRUE)
  }))
  predhist$median = unlist(sapply(predhists, function(x) {
    median(x, na.rm = TRUE)
  }))
  if (doPlot) {
    opar = par(no.readonly = TRUE)
    par(mai = mai)
    bars = barplot(predhist$mean, horiz = TRUE, col = barcol, 
                   las = 1, cex.names = cex.names, main = tit, ...)
    text(x = 0, y = bars, labels = sapply(predhists, function(x) paste(c("n = ", 
                                                                         length(x)), collapse = "")), pos = 4)
    par(opar)
  }
  return(predhist)
}


predHist7 <- function (nmin = 1, nmax, max_res = 500, ylim = c(0, 1), par.mai = c(1.02, 0.86, 0.82, 0.02), auc = FALSE, normalized = FALSE, xlab = "Top n Models", 
          ylab = "Frequency of Features", main = "Ranking of Features in Top Models", 
          cex.names = 0.6, doPlot = TRUE, plotDistributionAtNmax = FALSE, barcol="orange") 
{
  sampling <- if (nmax - nmin + 1 < max_res) 
    nmin:nmax
  else round(seq(nmin, nmax, (nmax - nmin + 1)/max_res))
  sampling[length(sampling)] <- nmax
  H <- t(sapply(sampling, function(y) {
    predHist(y, orderLike = NA, doPlot = FALSE)
  }))
  if (auc) {
    H = apply(H, 2, cumsum)
    H = H/1:nrow(H)
  }
  if (normalized) {
    H = H/rowSums(H)
  }
  opar = par(no.readonly = TRUE)
  par(mai = par.mai)
  if (is.na(ylim)) {
    ylim = range(H)
  }
  if (doPlot) {
    if (plotDistributionAtNmax) {
      barplot(100 * (sort(H[nmax, ], decreasing = FALSE)), 
              horiz = TRUE, col = barcol, las = 1, cex.names = cex.names, 
              xlab = xlab, main = main)
    }
    else {
      plot(c(0, 1.4 * nmax), ylim, type = "n", xlab = xlab, 
           ylab = ylab, main = main, cex.main = 1.7, cex.axis = 1.5, 
           cex.lab = 1.5)
      if (ind_sig <= nmax) 
        abline(v = ind_sig, col = "grey")
      abline(h = 0.5, col = "lightgrey")
      for (i in 1:dim(H)[2]) lines(sampling, H[, i], col = rainbow(dim(H)[2])[i])
      text(nmax, H[dim(H)[1], ], labels = colnames(H), 
           col = rainbow(dim(H)[2]), pos = 4, cex = 1)
    }
  }
  par(opar)
}



### ROCR f_r beste Marker
###
roc_one <- function(predNameList = list(), targetname, X=X, tit="ROC") {
  aucs <- list()
  roc_list <- list()
  j = 1
  for (predNames in predNameList) {
    indRowNA <- which(apply(is.na(d.f(X[,predNames])),1,sum) > 0)
    if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
    T <- as.factor(X[-indRowNA,targetname])
    Y <- d.f(d.f(X[-indRowNA,predNames]))
    names(Y) <- predNames
    fit <- glm(T ~ ., data=d.f(Y), family=binomial)
    pred <- prediction( fit$fitted.values, T)
    perf <- performance(pred, "sens", "spec")
    slot(perf, "x.values") <- lapply(slot(perf, "x.values"), function(x){1-x})
    slot(perf, "x.name") <- "1-Specifity"
    plot(perf, col=j, add=(j!=1))
    roc_list[[j]] <- perf
    aucs[[j]] <- slot(performance(pred, "auc"), "y.values")[[1]]
    j <- j+1
  }
  title(tit)
  legend("bottomright", unlist(lapply(predNameList, function(x){paste(x, collapse=", ")})), lty=1, col=1:length(predNameList), cex=0.5)
  for (i in 1:length(aucs)) {
    text(paste("AUC =",round(aucs[[i]], 3)), x=0.75, y=0.8-i/25, pos = 4, col=i)
  }
  return(roc_list)
}

roc_two <- function(predNameList = list(), targetname, X=X, tit="ROC", smoo = FALSE) {
  aucs <- list()
  roc_list <- list()
  j = 1
  for (predNames in predNameList) {
    indRowNA <- which(apply(is.na(d.f(X[,predNames])),1,sum) > 0)
    if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
    T <- as.factor(X[-indRowNA,targetname])
    Y <- d.f(d.f(X[-indRowNA,predNames]))
    names(Y) <- predNames
    fit <- glm(T ~ ., data=d.f(Y), family=binomial)
    ro <- tryCatch(pROC::roc(T, fit$fitted.values, smooth=smoo), error=function(x) {pROC::roc(T, fit$fitted.values, smooth=!smoo)})
    plot(ro, col=j, add=(j!=1))
    roc_list[[j]] <- ro
    aucs[[j]] <- ro$auc
    j <- j+1
  }
  title(tit, line=3)
  legend("bottomright", unlist(lapply(predNameList, function(x){paste(shortenNames(x, reallyshort=TRUE), collapse=", ")})), lty=1, col=1:length(predNameList), inset=0.05, cex=0.6)
  for (i in 1:length(aucs)) {
    text(paste("AUC =",round(aucs[[i]], 3)), x=0.2, y=0.8-i/25, pos = 4, col=i)
  }
  return(roc_list)
}

rocplus <- function(predNameList = list(), targetname, X=X, tit="ROC", compare="CA125", smoo = FALSE) {
  #  predNameListplus = list(compare)
  #  for (pn in predNameList) {
  #     predNameListplus <- c(predNameListplus, list(pn), list(c(compare, pn)))
  #  }
  #  roc_two(predNameList, targetname, X, tit)
  #  
  
  aucs <- list()
  
  # ROC compare
  predNames <- compare
  indRowNA <- which(apply(is.na(d.f(X[,predNames])),1,sum) > 0)
  if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
  T <- as.factor(X[-indRowNA,targetname])
  Y <- d.f(d.f(X[-indRowNA,predNames]))  
  names(Y) <- predNames
  fit <- glm(T ~ ., data=d.f(Y), family=binomial)
  ro <- tryCatch(pROC::roc(T, fit$fitted.values, smooth=smoo), error=function(x) {pROC::roc(T, fit$fitted.values, smooth=!smoo)})
  plot(ro, col=1, add=FALSE)
  aucs[[1]] <- ro$auc
  
  # ROC andere  
  j <- 2
  for (predNames in predNameList) {
    # ohne compare
    indRowNA <- which(apply(is.na(d.f(X[,predNames])),1,sum) > 0)
    if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
    T <- as.factor(X[-indRowNA,targetname])
    Y <- d.f(d.f(X[-indRowNA,predNames]))
    names(Y) <- predNames
    fit <- glm(T ~ ., data=d.f(Y), family=binomial)
    ro <- tryCatch(pROC::roc(T, fit$fitted.values, smooth=smoo), error=function(x) {pROC::roc(T, fit$fitted.values, smooth=!smoo)})
    plot(ro, col=j, lty=1, add=TRUE)
    aucs[[2*(j-1)]] <- ro$auc
    
    
    # mit compare
    predNames <- c(compare, predNames)
    indRowNA <- which(apply(is.na(d.f(X[,predNames])),1,sum) > 0)
    if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
    T <- as.factor(X[-indRowNA,targetname])
    Y <- d.f(d.f(X[-indRowNA,predNames]))
    names(Y) <- predNames
    fit <- glm(T ~ ., data=d.f(Y), family=binomial)
    ro <- tryCatch(pROC::roc(T, fit$fitted.values, smooth=smoo), error=function(x) {pROC::roc(T, fit$fitted.values, smooth=!smoo)})
    plot(ro, col=j, lty=2, add=TRUE)
    aucs[[2*j-1]] <- ro$auc
    
    j <- j+1
  }
  
  title(tit, line=3)
  legend("bottomright",
         c(compare, unlist(lapply(predNameList, function(x){c(paste(x, collapse=", "), paste(c(paste(x, collapse=", "), compare), collapse=", "))}))),
         lty=c(1, rep(1:2, times=length(predNameList))),
         col=c(1, rep(2:(1+length(predNameList)), each=2)),
         inset=0.05,
         cex=0.5)
  for (i in 1:length(aucs)) {
    text(paste("AUC =",round(aucs[[i]], 3)), x=0.25, y=0.7-i/25, pos = 4, col=ceiling((i+0.5)/2))
  }  
}


### AUCS f_r Marker
###
aucs <- function(predNameList = list(), targetname, X=X, smoo = FALSE) {
  sapply(predNameList, function(x){
    indRowNA <- which(apply(is.na(d.f(X[,x])),1,sum) > 0)
    if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
    T <- as.factor(X[-indRowNA,targetname])
    Y <- d.f(d.f(X[-indRowNA,x]))
    names(Y) <- x
    fit <- glm(T ~ ., data=d.f(Y), family=binomial)
    ro <- roc(T, fit$fitted.values, smooth=smoo)
    ro$auc
  })
}

### RMSD for two vectors
rmsd <- function(x, y, ...) {
  if (all(is.na(x)) || all(is.na(y))) return(NA)
  r = 0
  if (length(x)!=length(y)) stop("x and y must have equal length")
  r = sqrt(sum((x-y)^2, ...)/length(x))
  return (r)
}


### Correlation box for two vectors. n-fold bootstrapping
### validation_type can be "BS" (bootstrapped validation), "LOOCV" (leave one out CV),
### "k-fold" (k-fold-CV)
cor_box <- function(x, y, validation_type="BS", k=100, seedj=1, target=c(), patients=c(), folds=NA, errorfun="corerror", cormethod="pearson", observers=c(), ...) {
  set.seed(seedj)
  out <- c()
  if (validation_type=="BS") {
    ## FOR BOOTSTRAP VALIDATION
    for (i in 1:k) {
      if (length(patients) > 0) {
        if (nlevels(target)>0) {
          bsind <- c()
          for (nl in 1:nlevels(target)) {
            ind_tmp <- which(target==(levels(target)[nl]))
            bsind_patient <- sample(unique(patients[ind_tmp]), length(unique(patients[ind_tmp])), replace=TRUE)
            for (p in bsind_patient) {
              bsind <- c(bsind, which(patients==p))
            }
          }      
        } else {
          bsind <- c()
          bsind_patient <- sample(unique(patients), length(unique(patients)), replace=TRUE)
          for (p in bsind_patient) {
            bsind <- c(bsind, which(patients==p))
          }
        }
      } else {
        bsind <- sample(1:length(x), length(x), replace=TRUE)
        
        # for factors: draw samples out of every class
        if (nlevels(target)>0) {
          bsind <- c()
          for (nl in 1:nlevels(target)) {
            ind_tmp <- which(target==(levels(target)[nl]))
            bsind <- c(bsind, sample(ind_tmp, length(ind_tmp), replace=TRUE))
          }      
        }
      }
      oob <- c(1:length(x))[is.na(match(1:length(x), bsind))]
      if (errorfun=="corerror") {
        out <- c(out, cor(x[oob], y[oob], method=cormethod, use="pairwise.complete.obs", ...))
      } else {
        out <- c(out, rmsd(x[oob], y[oob], na.rm=TRUE, ...))
      }
    }
  } else if (validation_type=="k-fold") {
    if (any(is.na(folds))) {
      if (length(patients) > 0) {
        if (nlevels(target)>0) {
          folds <- rep(0, dim(X)[1])
          max_k <- k
          for (nl in 1:nlevels(target)) {
            ind_tmp <- which(target==(levels(target)[nl]))
            max_k <- min(max_k, length(unique(patients[ind_tmp])))
          }
          for (nl in 1:nlevels(target)) {
            ind_tmp <- which(target==(levels(target)[nl]))
            folds_patient <- sample(rep(1:max_k, round((length(unique(patients[ind_tmp])))/max_k)))[1:(length(unique(patients[ind_tmp])))]
            for (pat in unique(patients)) {
              
              folds[which(patients==pat)] <- folds_patient[which(unique(patients)==pat)]
            }
          }
        } else {
          folds <- rep(0, dim(X)[1])
          folds_patient <- sample(rep(1:k, round((length(unique(patients)))/k)))[1:(length(unique(patients)))]
          for (pat in unique(patients)) {
            folds[which(patients==pat)] <- folds_patient[which(unique(patients)==pat)]
          }
        }
      } else {
        folds <- sample(rep(1:k, round((length(x))/k)))[1:(length(x))]
        if (nlevels(target)>0) {
          max_k <- k
          for (nl in 1:nlevels(target)) {
            ind_tmp <- which(target==(levels(target)[nl]))
            max_k <- min(max_k, length(ind_tmp))
          }
          for (nl in 1:nlevels(target)) {
            ind_tmp <- which(target==(levels(target)[nl]))
            folds[ind_tmp] <- sample(rep(1:max_k, ceiling((length(ind_tmp))/max_k))[1:(length(ind_tmp))])
          }
        }
      }
    }
    
    for (i in unique(folds)) {
      oob <- which(folds==i)
      if (errorfun=="corerror") {
        out <- c(out, cor(x[oob], y[oob], method=cormethod, use="pairwise.complete.obs", ...))
      } else {
        out <- c(out, rmsd(x[oob], y[oob], na.rm=TRUE, ...))
      }
    }
  } else if (validation_type=="LOOCV") {
    if (length(patients) > 0) {
      return (cor_box(x, y, validation_type="k-fold", k, seedj, target, patients=c(), folds=patients, errorfun, ...))
    } else {
      return (cor_box(x, y, validation_type="k-fold", k, seedj, target, patients, folds=1:length(x), errorfun, ...))
    }
  } else if (validation_type=="LOPOCV") {
    
    folds <- as.numeric(observers)
    for (i in unique(folds)) {
      oob <- which(folds==i)
      if (errorfun=="corerror") {
        out <- c(out, cor(x[oob], y[oob], method=cormethod, use="pairwise.complete.obs", ...))
      } else {
        out <- c(out, rmsd(x[oob], y[oob], na.rm=TRUE, ...))
      }
    }
  }
  return (out)
}

classifier <- function(form, data, family, labels, doScale=FALSE, ...) {
  if (doScale) {
    predNames. <- colnames(data)
    data <- scale(data)
    mu <- attr(data, "scaled:center")
    std <- attr(data, "scaled:scale")
    data <- d.f(data)
    for (i in 1:dim(data)[2]) {
      if (all(is.na(data[,i]))) {
        data[,i] <- rep(0, dim(data)[1])
      }
    }
    colnames(data) <- predNames.
  }
  
  if (family=="tobit") {
    fit <- AER::tobit(as.numeric(as.character(labels))~., data=data)
  } else if (family=="nnls") { ### calculate non-negative coefficients, but keep the structure of a "glm"
    require(nnls)
    fit <- glm(formula=form, data=data, family="gaussian", ...)
    X_tmp <- model.matrix(form, data=data)
    nnlsfit <- nnls(X_tmp, labels)
    coeffis_names <- names(fit$coefficients)
    fit$coefficients <- coef(nnlsfit)
    names(fit$coefficients) <- coeffis_names
  } else {
    fit <- glm(formula=form, data=data, family=family, ...)
    
  }
  
  if (doScale) {
    attr(fit, "mu") <- mu
    attr(fit, "std") <- std
  }
  return (fit)
}


### Errorbox for Marker - one prediction
###
errorbox_compartment <- function (Y, T, targetname, testsamples, errorfun = "aucerror", 
                                  thirdtestname = NA, thirdtesttrainingsamples = NA, thirdtesttestsamples = NA, 
                                  doScale = FALSE) 
{
  fittedvalues <- list()
  origlabels <- list()
  if (is.na(thirdtestname)) {
    thirdtestname <- targetname
  }
  else {
    thirdtestname <- names(testsamples)[grep(thirdtestname, 
                                             names(testsamples))]
  }
  if (errorfun == "corerror" || errorfun == "rmsderror") {
    bsfit <- classifier(as.numeric(as.character(T)) ~ ., 
                        data = d.f(Y), family = cl_family, labels = T, doScale = doScale)
    if (!is.na(thirdtesttrainingsamples)) {
      thirdtestbsfit <- classifier(as.numeric(as.character(T)) ~ 
                                     ., data = d.f(thirdtesttrainingsamples), family = cl_family, 
                                   labels = T, doScale = doScale)
    }
  }
  else if (length(unique(T)) == 2) {
    bsfit <- classifier(T ~ ., data = d.f(Y), family = "binomial", 
                        labels = T, doScale = doScale)
  }
  else {
    bsfit <- randomForest(y = as.factor(T), x = d.f(Y), ntree = 100, 
                          do.trace = F, importance = F, family = "binomial")
  }
  if (doScale) {
    predNames. <- colnames(Y)
    mu <- attr(bsfit, "mu")
    std <- attr(bsfit, "std")
    for (x in predNames.) {
      testsamples[, x] <- (testsamples[, x] - mu[x])/std[x]
    }
  }
  testsA <<- testsamples
  bsfitA <<- bsfit
  if (errorfun == "aucerror") {
    if (any(is.na(coef(bsfit))) | any(coef(bsfit)[2:length(coef(bsfit))] <= 
                                      0)) {
      fittedvalues <- append(fittedvalues, list(rep(NA, 
                                                    dim(testsamples)[1])))
      origlabels <- append(origlabels, list(as.factor(testsamples[, 
                                                                  thirdtestname])))
    }
    else {
      fittedvalues <- append(fittedvalues, list(predict(bsfit, 
                                                        newdata = testsamples, type = "response")))
      origlabels <- append(origlabels, list(as.factor(testsamples[, 
                                                                  thirdtestname])))
    }
  }
  else if (errorfun == "corerror" || errorfun == "rmsderror") {
    if (any(is.na(coef(bsfit)))) {
      fittedvalues <- append(fittedvalues, list(rep(NA, 
                                                    dim(testsamples)[1])))
      origlabels <- append(origlabels, list(testsamples[, 
                                                        thirdtestname]))
    }
    else {
      fittedvalues <- append(fittedvalues, list(predict(bsfit, 
                                                        newdata = testsamples, type = "response")))
      if (!is.na(thirdtesttestsamples)) {
        origlabels <- append(origlabels, list(predict(thirdtestbsfit, 
                                                      newdata = thirdtesttestsamples, type = "response")))
      }
      else {
        origlabels <- append(origlabels, list(testsamples[, 
                                                          thirdtestname]))
      }
    }
  }
  else {
    if (any(is.na(coef(bsfit)))) {
      fittedvalues <- append(fittedvalues, 0)
    }
    else {
      result <- predict(bsfit, newdata = testsamples)
      if (length(unique(result)) != length(unique(T))) {
        result <- as.numeric(result >= 0)
      }
      result <- factor(result)
      levels(result) <- levels(testsamples[, targetname])
      cm <- table(result, testsamples[, thirdtestname])
      acc <- sum(diag(prop.table(cm)))
      sens <- c()
      specs <- c()
      for (nk in 1:length(unique(T))) {
        sens <- c(sens, cm[nk, nk]/sum(cm[, nk]))
        specs <- c(specs, sum(cm[-nk, -nk])/sum(cm[, 
                                                   -nk]))
      }
      bacc <- mean(sens)
      if (length(unique(T)) == 2) {
        sens <- sens[1]
        specs <- specs[1]
      }
      error <- acc
      fittedvalues <- append(fittedvalues, error)
    }
  }
  return(list(origlabels, fittedvalues))
}

validateErrorboxCompartment <- function(origlabels, fittedvalues, smoo, errorfun, cormethod, observer = NA) {
  if (errorfun == "aucerror") {
    markerError <- c()
    for (i in 1:length(origlabels)){
      if(length(unique(fittedvalues[[i]]))>1 && length(unique(origlabels[[i]]))>1) {
        mE <- tryCatch(pROC::roc(origlabels[[i]], fittedvalues[[i]], smooth=smoo)$auc, error = function(x) {pROC::roc(origlabels[[i]], fittedvalues[[i]], smooth=!smoo)$auc})
      } else if (validation_type=="LOOCV") {
        mE <- tryCatch(pROC::roc(unlist(origlabels), unlist(fittedvalues), smooth=smoo)$auc, error = function(x) {pROC::roc(unlist(origlabels), unlist(fittedvalues), smooth=!smoo)$auc})
      } else {
        mE <- NA
      }
      markerError <- c(markerError, mE)
    }
  } else if (errorfun == "corerror") {
    markerError <- c()
    p.values <- c()
    rocs <- c()
    for (i in 1:length(origlabels)){
      mE <- abs(cor(origlabels[[i]], fittedvalues[[i]], method=cormethod))
      p.value <- tryCatch(cor.test(origlabels[[i]], fittedvalues[[i]], method=cormethod)$p.value, error=function(x){NA})
      markerError <- c(markerError, mE)
      p.values <- c(p.values, p.value)
      rocs_ <- c()
      if (!is.na(observer)) {
        for (ob in unique(observer)) {
          roc_ = tryCatch(pROC::roc(origlabels[[i]][which(observer==ob)]>=3, fittedvalues[[i]][which(observer==ob)], smooth=smoo), error = function(x) {tryCatch(pROC::roc(origlabels[[i]][which(observer==ob)]>=3, fittedvalues[[i]][which(observer==ob)], smooth=!smoo), error = function(x) {NA})})
          rocs_ <- c(rocs_, list(roc_))
        }
      } else {
        rocs_ = NA
      }
      rocs <- c(rocs, list(rocs_))
    }
    attr(markerError, "p.values") <- p.values
    attr(markerError, "rocs") <- rocs
  } else if (errorfun == "rmsderror") {
    markerError <- c()
    for (i in 1:length(origlabels)){
      mE <- rmsd(origlabels[[i]], fittedvalues[[i]], na.rm=TRUE)
      markerError <- c(markerError, mE)
    }
  } else {
    markerError <- unlist(fittedvalues)
  }
  
  if (sum(is.na(markerError)) > 0.5*length(markerError)) {
    markerError <- rep(NA, length(markerError))
    p.values <- rep(NA, length(markerError))
    attr(markerError, "p.values") <- p.values
  }
  return(markerError)
}

### Errorbox for Marker
### 
# errorfun: "aucerror": optimization for AUC, "prederror": optimization for prediction
# validation_type can be "BS" (bootstrapped validation), "LOOCV" (leave one out CV, train on x-1 patients, calculate error on 1 patient),
# "w-fold" (w-fold-CV), "LOPOCV" (leave one patient out CV, train x times on x-1 patients and predict each time the 1 patient. In the end, calculate the error on all predictions).
errorbox <- function (predNames = c(), targetname, X, errorfun = "aucerror", 
                      seedj = 1, smoo = FALSE, randomtest = FALSE, validation_type = "BS", 
                      k = 100, thirdtestname = NA, thirdtestnewtrain = FALSE, bagging = FALSE, 
                      folds = NA, patientcol = -1, cormethod = "pearson", doScale = FALSE, 
                      observercol = 1) 
{
  if (length(predNames) < 1) 
    return(c())
  set.seed(seedj)
  B_tmp <- B
  if (randomtest) {
    X[, targetname] <- sample(X[, targetname])
    B_tmp[, targetname] <- sample(B_tmp[, targetname])
  }
  set.seed(seedj)
  fittedvalues <- list()
  origlabels <- list()
  test_obs <- list()
  if (validation_type == "BS") {
    bsfits <<- c()
    for (i in 1:k) {
      oob_patients <- c()
      if (patientcol > 0) {
        if (nlevels(X[, targetname]) > 0) {
          bsind <- c()
          for (nl in 1:nlevels(X[, targetname])) {
            ind_tmp <- which(X[, targetname] == (levels(X[, 
                                                          targetname])[nl]))
            bsind_patient <- sample(unique(X[ind_tmp, 
                                             patientcol]), length(unique(X[ind_tmp, 
                                                                           patientcol])), replace = TRUE)
            for (p in bsind_patient) {
              bsind <- c(bsind, which(X[, patientcol] == 
                                        p))
            }
          }
        }
        else {
          bsind <- c()
          bsind_patient <- sample(unique(X[, patientcol]), 
                                  length(unique(X[, patientcol])), replace = TRUE)
          for (p in (bsind_patient)) {
            bsind <- c(bsind, which(X[, patientcol] == 
                                      p))
          }
        }
        oob_patients <- unique(X[, patientcol])[-which(unique(X[, 
                                                                patientcol]) %in% bsind_patient)]
      }
      else {
        bsind <- sample(1:dim(X)[1], dim(X)[1], replace = TRUE)
        if (nlevels(X[, targetname]) > 0) {
          bsind <- c()
          for (nl in 1:nlevels(X[, targetname])) {
            ind_tmp <- which(X[, targetname] == (levels(X[, 
                                                          targetname])[nl]))
            bsind <- c(bsind, sample(ind_tmp, length(ind_tmp), 
                                     replace = TRUE))
          }
        }
      }
      if (FALSE & length(oob_patients) > 0) {
        oob <- which(B_tmp[, patientcol] %in% oob_patients)
        testsamples <- B_tmp[oob, ]
        test_obs_tmp <- B_tmp[oob, observercol]
      }
      else {
        oob <- c(1:dim(X)[1])[is.na(match(1:dim(X)[1], 
                                          bsind))]
        testsamples <- X[oob, ]
        test_obs_tmp <- X[oob, observercol]
      }
      if (thirdtestnewtrain) {
        thirdtesttestsamples = X[oob, ]
      }
      else {
        thirdtesttestsamples = NA
      }
      T <- as.factor(X[bsind, targetname])
      Y <- d.f(d.f(X[bsind, predNames]))
      names(Y) <- predNames
      thirdtesttrainingsamples <- NA
      if (thirdtestnewtrain) {
        thirdtesttrainingsamples <- d.f(d.f(X[, thirdtestname])[bsind, 
                                                                ])
        names(thirdtesttrainingsamples) <- thirdtestname
      }
      compartment <- errorbox_compartment(Y, T, targetname, 
                                          testsamples, errorfun, thirdtestname, thirdtesttrainingsamples, 
                                          thirdtesttestsamples, doScale)
      origlabels <- append(origlabels, compartment[[1]])
      fittedvalues <- append(fittedvalues, compartment[[2]])
      test_obs <- append(test_obs, test_obs_tmp)
      bsfits <<- c(bsfits, list(bsfitA))
    }
  }
  else if (validation_type == "k-fold") {
    if (any(is.na(folds))) {
      if (patientcol > 0) {
        if (nlevels(X[, targetname]) > 0) {
          folds <- rep(0, dim(X)[1])
          max_k <- k
          for (nl in 1:nlevels(X[, targetname])) {
            ind_tmp <- which(X[, targetname] == (levels(X[, 
                                                          targetname])[nl]))
            max_k <- min(max_k, length(unique(X[ind_tmp, 
                                                patientcol])))
          }
          for (nl in 1:nlevels(X[, targetname])) {
            ind_tmp <- which(X[, targetname] == (levels(X[, 
                                                          targetname])[nl]))
            folds_patient <- sample(rep(1:max_k, round((length(unique(X[ind_tmp, 
                                                                        patientcol])))/max_k)))[1:(length(unique(X[ind_tmp, 
                                                                                                                   patientcol])))]
            for (pat in unique(X[, patientcol])) {
              folds[which(X[, patientcol] == pat)] <- folds_patient[which(unique(X[, 
                                                                                   patientcol]) == pat)]
            }
          }
        }
        else {
          folds <- rep(0, dim(X)[1])
          folds_patient <- sample(rep(1:k, round((length(unique(X[, 
                                                                  patientcol])))/k)))[1:(length(unique(X[, 
                                                                                                         patientcol])))]
          for (pat in unique(X[, patientcol])) {
            folds[which(X[, patientcol] == pat)] <- folds_patient[which(unique(X[, 
                                                                                 patientcol]) == pat)]
          }
        }
      }
      else {
        folds <- sample(rep(1:k, round((dim(X)[1])/k)))[1:(dim(X)[1])]
        if (nlevels(X[, targetname]) > 0) {
          max_k <- k
          for (nl in 1:nlevels(X[, targetname])) {
            ind_tmp <- which(X[, targetname] == (levels(X[, 
                                                          targetname])[nl]))
            max_k <- min(max_k, length(ind_tmp))
          }
          for (nl in 1:nlevels(X[, targetname])) {
            ind_tmp <- which(X[, targetname] == (levels(X[, 
                                                          targetname])[nl]))
            folds[ind_tmp] <- sample(rep(1:max_k, ceiling((length(ind_tmp))/max_k))[1:(length(ind_tmp))])
          }
        }
      }
    }
    for (i in unique(folds)) {
      oob <- which(folds == i)
      testsamples <- X[oob, ]
      if (thirdtestnewtrain) {
        thirdtesttestsamples = X[oob, ]
      }
      else {
        thirdtesttestsamples = NA
      }
      bsind <- which(folds != i)
      T <- as.factor(X[bsind, targetname])
      Y <- d.f(X[bsind, predNames])
      names(Y) <- predNames
      thirdtesttrainingsamples <- NA
      if (thirdtestnewtrain) {
        thirdtesttrainingsamples <- d.f(d.f(X[, thirdtestname])[bsind, 
                                                                ])
        names(thirdtesttrainingsamples) <- thirdtestname
      }
      compartment <- errorbox_compartment(Y, T, targetname, 
                                          testsamples, errorfun, thirdtestname, thirdtesttrainingsamples, 
                                          thirdtesttestsamples, doScale)
      origlabels <- append(origlabels, compartment[[1]])
      fittedvalues <- append(fittedvalues, compartment[[2]])
    }
  }
  else if (validation_type == "LOOCV") {
    if (patientcol > 0) {
      return(errorbox(predNames, targetname, X, errorfun, 
                      seedj, smoo, randomtest, validation_type = "k-fold", 
                      k, thirdtestname, thirdtestnewtrain, bagging, 
                      folds = X[, patientcol], patientcol = -1, cormethod = cormethod))
    }
    else {
      return(errorbox(predNames, targetname, X, errorfun, 
                      seedj, smoo, randomtest, validation_type = "k-fold", 
                      k, thirdtestname, thirdtestnewtrain, bagging, 
                      folds = 1:dim(X)[1], patientcol, cormethod = cormethod))
    }
  }
  else if (validation_type == "LOPOCV") {
    X_tmp <- d.f(X[, predNames])
    names(X_tmp) <- predNames
    T_tmp <- X[, targetname]
    B$label_pred_cv <- rep(NA, dim(B)[1])
    for (pat in unique(X[, patientcol])) {
      T_all_fold = T_tmp[which(X[, patientcol] != pat)]
      X_all_fold = d.f(X_tmp[which(X[, patientcol] != pat), 
                             ])
      names(X_all_fold) <- predNames
      X_all_fold_train <- d.f(X_all_fold[, predNames])
      names(X_all_fold_train) <- predNames
      compartment <- errorbox_compartment(X_all_fold_train, 
                                          T_all_fold, targetname, B[which(B[, patientcol] == 
                                                                            pat), ], errorfun, NA, NA, NA, doScale)
      B$label_pred_cv[which(B[, patientcol] == pat)] <- unlist(compartment[[2]])
    }
    for (o in levels(B[, observercol])) {
      origlabels <- append(origlabels, list(B[which(B[, 
                                                      observercol] == o & B[, patientcol] %in% X[, 
                                                                                                 patientcol]), targetname]))
      fittedvalues <- append(fittedvalues, list(B$label_pred_cv[which(B[, 
                                                                        observercol] == o & B[, patientcol] %in% X[, 
                                                                                                                   patientcol])]))
    }
  }
  if (bagging) {
    out <- c()
    out$fittedvalues <- fittedvalues
    out$origvalues <- origlabels
    return(out)
  }
  markerError <- validateErrorboxCompartment(origlabels, fittedvalues, 
                                             smoo, errorfun, cormethod, observer = test_obs)
  attr(markerError, "CVmodels") = bsfits
  return(markerError)
}

### Errorbox for Marker for LOOCV and and correlation error:
### First, do the LOOCV prediction for all samples, Then, evaluate the correlations for 4 radiologists
### The errorbox has 4 values therefore.
# errorfun: only corerror allowed.
# validation_type: only "LOOCV" allowed
errorbox2 <- function(predNames = c(), targetname, X, errorfun="corerror", seedj = 1, randomtest = FALSE, validation_type="LOOCV", k=k, folds=NA, patientcol=-1, observercol=-1, cormethod="pearson") {
  if (length(predNames)<1) return(c())
  
  set.seed(seedj)
  if (randomtest) {
    X[,targetname] <- sample(X[,targetname])
  }
  
  set.seed(seedj) # fuer reproduzierbare Zufallszahlen
  
  X$orig_labels <- rep(NA, dim(X)[1])
  X$pred_labels <- rep(NA, dim(X)[1])
  
  if (validation_type=="k-fold") {
    ## FOR k-fold C-VALIDATION
    if (any(is.na(folds))) {
      if (patientcol > 0) {
        if (nlevels(X[,targetname])>0) {
          folds <- rep(0, dim(X)[1])
          max_k <- k
          for (nl in 1:nlevels(X[,targetname])) {
            ind_tmp <- which(X[,targetname]==(levels(X[,targetname])[nl]))
            max_k <- min(max_k, length(unique(X[ind_tmp,patientcol])))
          }
          for (nl in 1:nlevels(X[,targetname])) {
            ind_tmp <- which(X[,targetname]==(levels(X[,targetname])[nl]))
            folds_patient <- sample(rep(1:max_k, round((length(unique(X[ind_tmp,patientcol])))/max_k)))[1:(length(unique(X[ind_tmp,patientcol])))]
            for (pat in unique(X[,patientcol])) {
              folds[which(X[,patientcol]==pat)] <- folds_patient[which(unique(X[,patientcol])==pat)]
            }
          }
        } else {
          folds <- rep(0, dim(X)[1])
          folds_patient <- sample(rep(1:k, round((length(unique(X[,patientcol])))/k)))[1:(length(unique(X[,patientcol])))]
          for (pat in unique(X[,patientcol])) {
            folds[which(X[,patientcol]==pat)] <- folds_patient[which(unique(X[,patientcol])==pat)]
          }
        }
      } else {
        folds <- sample(rep(1:k, round((dim(X)[1])/k)))[1:(dim(X)[1])]
        if (nlevels(X[,targetname])>0) {
          max_k <- k
          for (nl in 1:nlevels(X[,targetname])) {
            ind_tmp <- which(X[,targetname]==(levels(X[,targetname])[nl]))
            max_k <- min(max_k, length(ind_tmp))
          }
          for (nl in 1:nlevels(X[,targetname])) {
            ind_tmp <- which(X[,targetname]==(levels(X[,targetname])[nl]))
            folds[ind_tmp] <- sample(rep(1:max_k, ceiling((length(ind_tmp))/max_k))[1:(length(ind_tmp))])
          }
        }
      }
    }
    
    for (i in unique(folds)) {
      
      #generiere Testset
      oob <- which(folds==i)
      testsamples <- X[oob,]
      
      # generiere Trainingsset
      bsind <- which(folds!=i)
      
      T <- as.factor(X[bsind,targetname])
      
      Y <- d.f(X[bsind,predNames])
      names(Y) <- predNames
      
      compartment <- errorbox_compartment(Y, T, targetname, testsamples, errorfun, doScale=doScale)
      
      X$orig_labels[which(folds==i)] <- compartment[[1]][[1]]
      X$pred_labels[which(folds==i)] <- compartment[[2]][[1]]
    }
  } else if (validation_type=="LOOCV") {
    if (patientcol>0) {
      return (errorbox2(predNames, targetname, X, errorfun, seedj, randomtest, validation_type="k-fold", k=k, folds=X[,patientcol], patientcol=-1, observercol=observercol, cormethod=cormethod))
    } else {
      return (errorbox2(predNames, targetname, X, errorfun, seedj, randomtest, validation_type="k-fold", k=k, folds=1:dim(X)[1], patientcol=patientcol, observercol=observercol, cormethod=cormethod))
    }
  }
  
  
  markerError <- c()
  if (observercol>0) {
    observers <- X[,observercol]
  } else {
    observers <- rep(1, dim(X)[1])
  }
  
  if (errorfun == "corerror") {
    p.values <- c()
    for (obs in levels(as.factor(observers))) {
      # local CDEIS correlation
      a <- X$orig_labels[which(observers==obs)]
      b <- X$pred_labels[which(observers==obs)]
      mE <- abs(cor(a, b, method=cormethod))
      p.value <- tryCatch(cor.test(a, b, method=cormethod)$p.value, error=function(x){NA})
      markerError <- c(markerError, mE)
      p.values <- c(p.values, p.value)
    }
    attr(markerError, "p.values") <- p.values
  } else if (errorfun == "rmsderror") {
    for (obs in levels(as.factor(observers))) {
      # local CDEIS correlation
      a <- X$orig_labels[which(observers==obs)]
      b <- X$pred_labels[which(observers==obs)]
      mE <- rmsd(a, b, na.rm=TRUE)
      markerError <- c(markerError, mE)
    }
  } else {
    markerError <- unlist(X$pred_labels)
  }    
  return(markerError)
}


### AUC bootstrapping for one marker
### 
aucbootstrap <- function(predNames = c(), targetname, X, seedj = 1, smoo = FALSE, disp=TRUE, randomtest = FALSE, folds = NA) {
  if (length(predNames)<1) return(c())
  
  set.seed(seedj)
  if (randomtest) {
    X[,targetname] <- sample(X[,targetname])
  }
  
  if (is.na(folds)) {
    set.seed(seedj)
  } else {
    k <- length(unique(folds))
  }
  
  perfs <- list();
  m_auc <- c()
  for (i in 1:k) {
    if (is.na(folds)) {
      if (nlevels(X[,targetname])>0) {
        bsind <- c()
        for (nl in 1:nlevels(X[,targetname])) {
          ind_tmp <- which(X[,targetname]==(levels(X[,targetname])[nl]))
          bsind <- c(bsind, sample(allind, length(allind), replace=TRUE))
        }      
      } else {
        bsind <- sample(allind, length(allind), replace=TRUE)
      }
      oob <- allind[is.na(match(allind, bsind))]
    } else {
      oob <- which(folds==i)
      bsind <- which(folds!=i)
    }
    testsamples <- X[oob,]
    T <- as.factor(X[bsind,targetname])
    
    Y <- d.f(d.f(X[bsind,predNames]))
    names(Y) <- predNames
    bsfit <- glm(T ~ ., data=d.f(Y), family=binomial)
    
    fittedvalues <- predict(bsfit, newdata=testsamples, type="response")
    origlabels <- as.factor(testsamples[,targetname])
    
    ro <- tryCatch(roc(origlabels, fittedvalues, smooth=smoo), error=function(x) {roc(origlabels, fittedvalues, smooth=!smoo)})
    if (disp) {
      plot(ro, col="gray", add=(i!=1))
    }
    
    perfs <- c(perfs, list(ro))
    
    m_auc <- c(m_auc, ro$auc)
  }
  perfs <- avg.roc(perfs)
  if (disp) {
    plot(perfs, col="red", add=TRUE)
    
    #indRowNA <- which(apply(is.na(d.f(X[,predNames])),1,sum) > 0)
    #  if (length(indRowNA) == 0 ) indRowNA = dim(X)[1] + 1
    #  T <- as.factor(X[-indRowNA,targetname])
    #  Y <- d.f(d.f(X[-indRowNA,predNames]))
    #  names(Y) <- predNames
    #  fit <- glm(T ~ ., data=d.f(Y), family=binomial)
    #  pred <- prediction( fit$fitted.values, T)
    #  perf <- performance(pred, "sens", "spec")
    #  slot(perf, "x.values") <- lapply(slot(perf, "x.values"), function(x){1-x})
    #  slot(perf, "x.name") <- "1-Specifity"
    #  plot(perf, col="red", add=TRUE)
    
    title(paste(c(as.character(k), "-fold bootstrap AUC of a marker"), collapse=""), line=3)
    legend("bottomright",
           c(paste(shortenNames(predNames, reallyshort=TRUE), collapse=", "), "mean"),
           lty=c(1),
           col=c("gray", "red"),
           inset=0.05,
           cex=0.6)
    text(paste("Median AUC =",round(median(m_auc, na.rm=TRUE), 3)), x=0.3, y=0.6-1/25, pos = 4, col="red")
    text(paste("Mean AUC =",round(mean(m_auc, na.rm=TRUE), 3)), x=0.3, y=0.6-2/25, pos = 4, col="black")
    text(paste("Spec at 75% Sens =",round(100*Sens4Spec(perfs, 0.75), 1), "%"), x=0.3, y=0.6-3/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at 90% Sens =",round(100*Sens4Spec(perfs, 0.9), 1), "%"), x=0.3, y=0.6-4/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at 98% Sens =",round(100*Sens4Spec(perfs, 0.98), 1), "%"), x=0.3, y=0.6-5/25, pos = 4, cex=0.6, col="blue")
  }  
  return(m_auc)
  #return(perfs)
}


##
## Sensitivity at specific Specificity for a given performance
## (What is the Sens. for 0.75 Spec in the AUC?)
Sens4Spec <- function(perf, se) {
  #sp <- slot(perf, "y.values")[[1]][tail(which(slot(perf, "x.values")[[1]]<=(1-se)),1)]
  sp <- perf$specificities[rev(which(perf$sensitivities>=(se)))[1]]
  return (sp)
}

###
### 
###
delete_successful_lsf <- function(foldername=".") {
  lsf_files <- dir(foldername, pattern='lsf.*');
  
  if (length(lsf_files) > 0) {
    for (i in 1:length(lsf_files)) {
      words = read.table(lsf_files[i], sep="\n", stringsAsFactors=FALSE, nrows=30)
      if (any(words == "Successfully completed.")) {
        unlink(lsf_files[i]);
      }
    }
  }
}

###
### 
###
resubmit_unsuccessful_lsf <- function(foldername, bsub_cmd = "bsub -W 08:00 ") {
  lsf_files <- dir(foldername, pattern='lsf.*');
  
  if (length(lsf_files) > 0) {
    for (i in 1:length(lsf_files)) {
      words = read.table(lsf_files[i], sep="\n", stringsAsFactors=FALSE, nrows=30)
      if (!any(words == "Successfully completed.")) {
        ind <- which(words == "Your job looked like:")
        cmd <- paste(c(bsub_cmd, " ", words[ind+2,]), collapse="")
        system(cmd, intern = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
        unlink(lsf_files[i]);
      }
    }
  }
}


###
### Brute Force Marker Screening Index erstellung
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

CreateIndices <- function(n, p) {
  out <- list()
  for (i in 1:p) {
    out <- append(out, CreateIndices_fast(n, i))
  }
  return (out)
}


###
### Brute Force Marker Screening Index erstellung
###
GetIndex <- function(n, p, indices) {
  if (length(indices)==1) {
    return (indices)
  } else {
    out <- CreateIndices_rek(n, p-1)
    ind <- length(out)
    sizes <- sapply(out, length)
    max_size <- max(sizes)
    ind_expand <- which(sizes==max_size)
    for (i in 1:n) {
      for (j in ind_expand) {
        if (i<head(out[[j]],1)) {
          ind <- ind + 1
          out[[ind]] <- c(i, out[[j]])
        }
      }
    }
    return (out)
  }
}


### Apply a classifier
### 
# applyClassifier: "aucerror": optimization for AUC, "prederror": optimization for prediction
applyClassifier <- function(X_train, Y_train, X_test, Y_test, X_add = c(), Y_add = c(), errorfun="aucerror" , classnames = c(), seedj = 1, smoo=FALSE, cutoff=0.5, cormethod="pearson") {
  set.seed(seedj) # fuer reproduzierbare Zufallszahlen
  
  fittedvalues <- c()
  origlabels <- c()
  
  
  
  if (errorfun == "corerror") {
    bsfit <- glm(as.numeric(as.character(Y_test)) ~ ., data=d.f(X_train), family=gaussian)
  } else if (length(unique(Y_train))==2) {
    bsfit <- glm(Y_train ~ ., data=X_train, family=binomial)
  } else {
    bsfit <- randomForest(y=as.factor(Y_train), x=X_train, ntree=100, do.trace=F, importance=F, family="binomial")
  }
  
  if (errorfun == "aucerror") {
    # speichere die fitted values und labels der OOB samples
    fittedvalues <- predict(bsfit, newdata=X_test, type="response")
    origlabels <- Y_test
    #} else if (errorfun == "corerror") {
    #     bsfit <- glm(as.numeric(as.character(Y_test)) ~ ., data=d.f(X_train), family=gaussian)
    #     if (any(is.na(coef(bsfit)))) {
    #       fittedvalues <- rep(NA, dim(Y_test)[1])
    #    	 origlabels <- Y_test
    #   	 } else {
    #   	   fittedvalues <- predict(bsfit, newdata=X_test, type="response")
    #       origlabels <- Y_test
    #     }
  } else {
    # klassifiziere und speichere den error
    if (length(unique(Y_train))==2 || errorfun == "corerror") {
      pre_result <- predict(bsfit, newdata=X_test, type="response")
      result <- factor(pre_result >= cutoff)
      result <- as.numeric(result)-1
    } else {
      pre_result <- apply(predict(bsfit, newdata=X_test, type="prob"), 1, max)
      result <- predict(bsfit, newdata=X_test, type="response")
    }
    
    #cm = table(Y_test, result) 
    ## absolute error
    #if (dim(cm)[2]==1) {
    #  error <- cm[1,1]/sum(cm)
    #} else {
    #  error <- (cm[1,2]+cm[2,1])/sum(cm)
    #}
    ## balanced error
    #if (dim(cm)[2]==1) {
    #  error <- 0.5
    #} else {
    #  error <- 1-mean(cm[1,2]/(cm[1,2]+cm[1,1]), cm[2,1]/(cm[2,1]+cm[2,2]))
    #}
    
    #fittedvalues <- error
  }
  
  
  
  if (errorfun == "aucerror") {
    ro <- roc(origlabels, fittedvalues, smooth=smoo)
    markerError <- ro$auc
    
    plot(ro, col=1, lty=1)
    title("Validation of classifier", line=3, sub=paste(c("Validation samples (", as.character(length(which(Y_test==levels(Y_test)[1]))),"|", as.character(length(which(Y_test==levels(Y_test)[2]))), ") have never been seen in analysis"), collapse=""))
    legend("bottomright",
           c(paste(names(X_train), collapse=", ")),
           lty=c(1),
           inset=0.05,
           cex=0.5)
    auc = ro$auc
    
    text(paste("AUC =",round(auc, 3)), x=0.25, y=0.7-1/25, pos = 4, col=ceiling((1+0.5)/2))
    text(paste("Spec at 75% Sens =",round(100*Sens4Spec(ro, 0.75), 1)), x=0.25, y=0.6-2/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at 90% Sens =",round(100*Sens4Spec(ro, 0.9), 1)), x=0.25, y=0.6-3/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at 98% Sens =",round(100*Sens4Spec(ro, 0.98), 1)), x=0.25, y=0.6-4/25, pos = 4, cex=0.6, col="blue")
    show(ro$sensitivities)
    show(ro$specificities)
    return(markerError)
    #} else if (errorfun == "corerror") {
    #      markerError <- abs(cor(origlabels, fittedvalues, method=cormethod))
    #      return(markerError)
    #    } else if (errorfun == "rmsderror") {
    #      markerError <- rmsd(origlabels[[i]], fittedvalues[[i]], na.rm=TRUE)
    #      return (markerError)
  } else {
    out <- list()
    out <- append(out, list(pre_result))
    out <- append(out, list(result))
    
    Y_test <- factor(Y_test)
    levels(Y_test) <- sapply(levels(Y_test), function(x){substr(x, 1,13)})
    result <- factor(result)
    levels(result) <- sapply(levels(Y_test), function(x){substr(x, 1,13)})
    cm <- table(result, Y_test, dnn=c("pred", "truth")) # truth top (cols), pred left (rows)
    tmpdimnamesnames <- names(dimnames(cm))
    if (length(classnames)==length(dimnames(cm)[[1]])) {
      tmpdimnames <- classnames
      dimnames(cm) <- list(tmpdimnames, tmpdimnames)
    } else {
      tmpdimnames <- dimnames(cm)[[1]]
    }
    
    # Acc
    acc <- sum(diag(cm))/sum(cm)
    
    # Sens, Spec
    sens <- c()
    specs <- c()
    for (nk in 1:nlevels(Y_test)) {
      sens <- c(sens, cm[nk, nk]/sum(cm[,nk]))
      specs <- c(specs, sum(cm[-nk,-nk])/sum(cm[,-nk]))
    }
    
    # balanced Acc ITS NOT CLEAR, WHAT BALANCED ACCURACY MEANS FOR MORE THAN TWO CLASSES
    bacc <- mean(sens)
    
    # balanced Acc ITS NOT CLEAR, WHAT BALANCED ACCURACY MEANS FOR MORE THAN TWO CLASSES
    #bacc <- mean((sens+specs)/2)
    
    
    cm <- cbind(rbind(cm, "total"=colSums(cm)), "total"=c(rowSums(cm), sum(cm)))
    cm <- rbind(cm, "sensitivity [%]"=c(round(100*sens), NA))
    cm <- rbind(cm, "specificity [%]"=c(round(100*specs), NA))
    
    ## classify additional samples
    if (length(X_add)>0) {
      if (length(unique(Y_train))==2) {
        pre_result <- predict(bsfit, newdata=X_add, type="response")
        result <- factor(pre_result >= cutoff)
        result <- as.numeric(result)-1
      } else {
        pre_result <- apply(predict(bsfit, newdata=X_add, type="prob"), 1, max)
        result <- predict(bsfit, newdata=X_test, type="response")
      }
      out_add <- matrix(rep(0,length(tmpdimnames)*length(unique(T_add_raw))), ncol=length(unique(T_add_raw)), dimnames=list(c(), unique(T_add_raw)))
      for (i in 1:dim(X_add)[1]) {
        out_add[result[i]+1, which(dimnames(out_add)[[2]]==T_add_raw[i])] <- out_add[result[i]+1, which(dimnames(out_add)[[2]]==T_add_raw[i])] + 1
      }       
      cm <- cbind(cm, rbind(out_add, matrix(rep(NA, 3*length(unique(T_add_raw))), ncol=length(unique(T_add_raw)))))
    }
    
    names(dimnames(cm)) = tmpdimnamesnames
    
    textplot(cm, mar=c(0,0,0,0), cmar=1, cex=0.5)
    title("Validation Set Prediction, Cols: truth, Rows: pred", sub=paste(c("cutoff = ", as.character(cutoff), " Accuracy = ", as.character(round(1000*acc)/10), "%", ", balanced Accuracy = ", as.character(round(1000*bacc)/10), "%"), collapse=""))
    #markerError <- unlist(fittedvalues)
    return(out)
  }          
}

shortenNames <- function(namelist, reallyshort=FALSE) {
  if (reallyshort) {
    return (sapply(namelist, function(x){strsplit(x, "_")[[1]][1]}, USE.NAMES=FALSE))
  } else {
    return (sapply(namelist, function(x){l = strsplit(x, "_")[[1]];  if (length(l)>1) {paste(c(l[1], substr(l[2],1,1)), collapse="_")} else {l[1]}}, USE.NAMES=FALSE))
  }
}

myImagePlot <- function(x, ...){
  p.old <- par(no.readonly=TRUE)
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  cexaxis <- 0.7
  ColorRamp <- NA
  pmar <- c(5,5,2.5,2)
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
    if( !is.null(Lst$cex.axis) ){
      cexaxis <- Lst$cex.axis
    }
    if( !is.null(Lst$col) ){
      ColorRamp <- Lst$col
    }
    if( !is.null(Lst$mar) ){
      pmar <- Lst$mar
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  if (any(is.na(ColorRamp))) {
    ColorRamp <- rgb( seq(0,1,length=256),  # Red
                      seq(0,1,length=256),  # Green
                      seq(1,0,length=256))  # Blue
  }
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = pmar)
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=cexaxis, las=2)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=cexaxis)
  
  # Color Scale
  par(mar = c(pmar[1], 2.5, pmar[3:4]))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
  par(p.old)
}

# bagging of classifiers: Take a set of models (list of predName vectors), train
# each classifier on X_train and Y_train and predict X_test. Use Majority vote
# for the final prediction. Compare it with Y_test, if available.
baggingClassifiers <- function(predNameList, X_train, Y_train, X_test, Y_test, X_add = c(), Y_add = c(), errorfun="aucerror" , classnames = c(), seedj = 1, smoo=FALSE, cutoff=0.5, cormethod="pearson") {
  bag <- c()
  for (i in 1:length(predNameList)) {
    predNames_ <- predNameList[[i]]
    X_train_ <- as.data.frame(X_train[,predNames_]); names(X_train_) <- predNames_
    X_test_ <- as.data.frame(X_test[,predNames_]); names(X_test_) <- predNames_
    out <- applyClassifier(X_train_, Y_train, X_test_, Y_test, X_add, Y_add, errorfun, classnames, seedj, smoo, cutoff)
    bag <- rbind(bag, (as.numeric(as.character(out[[1]]))))
  }
  prediction <- apply(bag, 2, mean)
  
  if (errorfun == "aucerror") {
    markerError <- tryCatch(pROC::roc(Y_test, prediction, smooth=smoo)$auc, error = function(x) {pROC::roc(Y_test, prediction, smooth=!smoo)$auc})
  } else if (errorfun == "corerror") {
    markerError <- abs(cor(Y_test, prediction, method=cormethod))
  } else if (errorfun == "rmsderror") {
    markerError <- rmsd(Y_test, prediction, na.rm=TRUE)
  } else {
    markerError <- unlist(prediction)
  }    
  
  out <- c(list(markerError), list(Y_train), list(prediction))
  return(out)  
}

errorbox_bagging <- function(predNameList, targetname, X=X, errorfun="aucerror" , seedj = 1, smoo = FALSE, randomtest = FALSE, validation_type="BS", k=100, cormethod="pearson") {
  fittedvalues <- list()
  origlabels <- list()
  for (i in 1:length(predNameList)) {
    predNames_ <- predNameList[[i]]
    out <- errorbox(predNames_, targetname, X, errorfun, seedj, smoo, randomtest, validation_type, k, bagging=TRUE)
    if (i==1) {
      fittedvalues <- out$fittedvalues
      origlabels <- out$origvalues
    } else {
      fittedvalues <- lapply(1:length(fittedvalues), function(x) {fittedvalues[[x]]+out$fittedvalues[[x]]})
      origlabels <- lapply(1:length(origlabels), function(x) {origlabels[[x]]+out$origvalues[[x]]})
    }
  }
  
  fittedvalues <- lapply(fittedvalues, function(x) {x/length(predNameList)})
  origlabels <- lapply(origlabels, function(x) {x/length(predNameList)})
  
  if (errorfun == "aucerror") {
    markerError <- c()
    for (i in 1:length(origlabels)){
      if(length(unique(fittedvalues[[i]]))>1 && length(unique(origlabels[[i]]))>1) {
        mE <- tryCatch(pROC::roc(origlabels[[i]], fittedvalues[[i]], smooth=smoo)$auc, error = function(x) {pROC::roc(origlabels[[i]], fittedvalues[[i]], smooth=!smoo)$auc})
      } else {
        mE <- NA
      }
      markerError <- c(markerError, mE)
    }
  } else if (errorfun == "corerror") {
    markerError <- c()
    for (i in 1:length(origlabels)){
      mE <- abs(cor(origlabels[[i]], fittedvalues[[i]], method=cormethod))
      markerError <- c(markerError, mE)
    }
  } else if (errorfun == "rmsderror") {
    markerError <- c()
    for (i in 1:length(origlabels)){
      mE <- rmsd(origlabels[[i]], fittedvalues[[i]], na.rm=TRUE)
      markerError <- c(markerError, mE)
    }
  } else {
    markerError <- unlist(fittedvalues)
  }    
  return(markerError)
}

eval_distribution <- function(numeric_var, threshold, factor_var) {
  e <- c()
  m <- contingency(numeric_var, threshold, factor_var)
  Accur <- (m[1,1] + m[2,2]) / (sum(m) + 0.0000000001)
  Spec <- m[2,2] / (m[2,2] + m[2,1] + 0.0000000001)
  TPrate <- m[1,1] / (m[1,1] + m[2,1] + 0.0000000001)
  TNrate <- m[2,2] / (m[2,2] + m[1,2] + 0.0000000001)
  fisher <- fisher.test(m)
  e <- data.frame(Accur, Spec, TPrate, TNrate, P.Value = fisher$p.value)
  return (e)
}

eval_distributions <- function(numeric_var, factor_var, withplot = TRUE) {
  n1 <- names(numeric_var)
  n2 <- names(factor_var)
  factor_var <- as.factor(factor_var[[1]])
  e <- c()
  mi <- min(numeric_var, na.rm=T)
  ma <- max(numeric_var, na.rm=T)
  thresholds <- seq(mi, ma, 0.5)
  for (i in thresholds) {
    e <- rbind(e, eval_distribution(numeric_var, i, factor_var))
  }
  if (withplot) {
    tit = paste(c("Threshold in ", n1, " fits ", n2), collapse="")
    subt = "Plot calc. from contingency tables with given threshold. P: Fisher Exact Test."
    xlabel = "Threshold"
    plot(thresholds, e[,1], type="l", ylim=c(0,1), main=tit, sub=subt, xlab=xlabel, ylab="1: good, 0:bad (expt P-Value)", cex.main=1, cex.sub=0.5)
    for (i in 2:dim(e)[2]) {
      lines(thresholds, e[,i], col=i)
    }
    legend("topright", names(e), lty=1, col=1:dim(e)[2])
  }
  return (e)  
}

contingency <- function(numeric_var, threshold, factor_var) {
  f <- levels(factor_var)
  m1 <- matrix(c(0,0,0,0), ncol=2)
  m1[1,1] <- sum(factor_var[which(numeric_var>=threshold)]==f[2], na.rm=TRUE)
  m1[1,2] <- sum(factor_var[which(numeric_var>=threshold)]==f[1], na.rm=TRUE)
  m1[2,1] <- sum(factor_var[which(numeric_var<threshold)]==f[2], na.rm=TRUE)
  m1[2,2] <- sum(factor_var[which(numeric_var<threshold)]==f[1], na.rm=TRUE)
  
  m2 <- matrix(c(0,0,0,0), ncol=2)
  m2[1,1] <- sum(factor_var[which(numeric_var>=threshold)]==f[1], na.rm=TRUE)
  m2[1,2] <- sum(factor_var[which(numeric_var>=threshold)]==f[2], na.rm=TRUE)
  m2[2,1] <- sum(factor_var[which(numeric_var<threshold)]==f[1], na.rm=TRUE)
  m2[2,2] <- sum(factor_var[which(numeric_var<threshold)]==f[2], na.rm=TRUE)
  
  Accur1 <- (m1[1,1] + m1[2,2]) / (sum(m1) + 0.0000000001)
  Accur2 <- (m2[1,1] + m2[2,2]) / (sum(m2) + 0.0000000001)
  
  if (Accur1 > Accur2) return (m1)
  return (m2)  
}

validationAUC <- function(X_train, Y_train, X_test, Y_test, cutoff=-1, crossvalidated=TRUE, subt="", k=k, folds=NA) {
  if (cutoff<0) {
    findCutoff <- TRUE
  } else {
    findCutoff <- FALSE
  }
  
  cutoffs <- (0:10)/10
  values <- c()
  values$accs <- list()
  values$specs <- list()
  values$senss <- list()
  for (i in 1:length(cutoffs)) {
    values$accs[[i]] <- list()
    values$specs[[i]] <- list()
    values$senss[[i]] <- list()
  }
  if (crossvalidated) {
    if (is.na(folds)) {
      set.seed(1)
    } else {
      k <- length(unique(folds))
    }
    perfs <- list()
    m_auc <- c()
    accs <- c()
    for (i in 1:k) {
      if (is.na(folds)) {
        if (nlevels(Y_test)>0) {
          bsind <- c()
          for (nl in 1:nlevels(Y_test)) {
            ind_tmp <- which(Y_test==(levels(Y_test)[nl]))
            bsind <- c(bsind, sample(ind_tmp, length(ind_tmp), replace=TRUE))
          }      
        } else {
          bsind <- sample(1:length(Y_test), length(Y_test), replace=TRUE)
        }
        oob <- c(1:length(Y_test))[is.na(match(1:length(Y_test), bsind))]
      } else {
        oob <- which(folds==i)
        bsind <- which(folds!=i)
      }
      
      X_train_ <- data.frame(X_train[bsind,])
      names(X_train_) <- names(X_train)
      bsfit <- glm(Y_train[bsind] ~ ., data=X_train_, family=binomial)
      probabilities <- predict(bsfit, newdata=X_test, type="response")
      probabilities_cutoff05 <- as.factor(as.numeric(probabilities>=cutoff))
      
      for (c in 1:length(cutoffs)) {
        probabilities_cutoff <- as.factor(as.numeric(probabilities>=cutoffs[c]))
        values$accs[[c]] <- append(values$accs[[c]], length(which((Y_test==0 & probabilities_cutoff==0)|(Y_test==1 & probabilities_cutoff==1)))/length(Y_test))
        values$specs[[c]] <- append(values$specs[[c]], length(which(Y_test==0 & probabilities_cutoff==0))/length(Y_test[which(Y_test==0)]))
        values$senss[[c]] <- append(values$senss[[c]], length(which(Y_test==1 & probabilities_cutoff==1))/length(Y_test[which(Y_test==1)]))
      }
      
      accs <- c(accs, length(which((Y_test==0 & probabilities_cutoff05==0)|(Y_test==1 & probabilities_cutoff05==1)))/length(Y_test))
      ro <- roc(Y_test, probabilities, smooth=FALSE)
      plot(ro, col="gray", add=(i!=1))
      perfs <- c(perfs, list(ro))
      m_auc <- c(m_auc, ro$auc)
    }
    perfs <- avg.roc(perfs)
    plot(perfs, col="red", add=TRUE, sub=subt)
    title(paste(c("Bootstrapped AUC of model ", paste(names(X_train), collapse=", ")), collapse=""), line=3)
    text(paste("Median AUC =",round(median(m_auc, na.rm=TRUE), 3)), x=0.3, y=0.6-1/25, pos = 4, cex=0.6, col="red")
    text(paste("Mean AUC =",round(mean(m_auc, na.rm=TRUE), 3), " +-", round(sqrt(var(m_auc, na.rm=TRUE)), 3)), x=0.3, y=0.6-2/25, pos = 4, cex=0.6, col="black")
    text(paste("Spec at 75% Sens =",round(100*Sens4Spec(perfs, 0.75), 1), "%"), x=0.3, y=0.6-3/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at 90% Sens =",round(100*Sens4Spec(perfs,  0.9), 1), "%"), x=0.3, y=0.6-4/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at 98% Sens =",round(100*Sens4Spec(perfs, 0.98), 1), "%"), x=0.3, y=0.6-5/25, pos = 4, cex=0.6, col="blue")
    text(paste("Acc at cutoff ", as.character(round(1000*cutoff[ind])/1000), " = ", round(100*mean(accs, na.rm=TRUE), 1), " +-", round(100*sqrt(var(accs, na.rm=TRUE)), 1), "%"), x=0.3, y=0.6-6/25, pos = 4, cex=0.6, col="blue")
    text(subt, x=0.5, y=0, pos=3)
    
    accus <- sapply(1:length(cutoffs), function(x){mean(unlist(values$accs[[x]]), na.rm=TRUE)})
    specs <- sapply(1:length(cutoffs), function(x){mean(unlist(values$specs[[x]]), na.rm=TRUE)})
    senss <- sapply(1:length(cutoffs), function(x){mean(unlist(values$senss[[x]]), na.rm=TRUE)})
    
    #lines(cutoffs, accs, col="green", lty=2)
    #lines(cutoffs, specs, col="blue", lty=3)
    #lines(cutoffs, senss, col="purple", lty=4)
    #legend(0.4, 0.2, c("Acc", "Spec", "Sens"), col=c("green", "blue", "purple"), lty=c(2,3,4), cex=0.6)
    
    # draw the point on the AUC with the best Accuracy 
    ind <- which.max(accus)
    text("+", x=specs[ind], y=senss[ind], pos = 4, cex=1.5, col="red")
    text(paste(c("(Acc = ", as.character(round(100*accus[ind])/100), ")"), collapse=""), x=specs[ind], y=senss[ind]-0.05, pos = 4, col="red")
    
    
  } else {
    bsfit <- glm(Y_train ~ ., data=X_train, family=binomial)
    probabilities <- predict(bsfit, newdata=X_test, type="response")
    
    for (c in 1:length(cutoffs)) {
      probabilities_cutoff <- as.factor(as.numeric(probabilities>=cutoffs[c]))
      values$accs[[c]] <- append(values$accs[[c]], length(which((Y_test==0 & probabilities_cutoff==0)|(Y_test==1 & probabilities_cutoff==1)))/length(Y_test))
      values$specs[[c]] <- append(values$specs[[c]], length(which(Y_test==0 & probabilities_cutoff==0))/length(Y_test[which(Y_test==0)]))
      values$senss[[c]] <- append(values$senss[[c]], length(which(Y_test==1 & probabilities_cutoff==1))/length(Y_test[which(Y_test==1)]))
    }
    ro <- roc(Y_test, probabilities, smooth=FALSE)
    plot(ro)
    title(paste(c("AUC of model ", paste(names(coef(bsfit))[-grep("\\(", names(coef(bsfit)))], collapse=", ")), collapse=""), line=3)
    
    prevalence <- 1-length(which(ro$response==ro$levels[1]))/length(ro$response)
    accus <- ro$specificities*(1-prevalence) + ro$sensitivities*(prevalence)
    #lines(ro$thresholds, ro$specificities*(1-prevalence) + ro$sensitivities*(prevalence), col="green", lty=2)
    #lines(ro$thresholds, ro$specificities, col="blue", lty=3)
    #lines(ro$thresholds, ro$sensitivities, col="purple", lty=4)
    #legend(0.4, 0.2, c("Acc", "Spec", "Sens"), col=c("green", "blue", "purple"), lty=c(2,3,4), cex=0.6)
    
    if(findCutoff) {
      # draw the point on the AUC with the best Accuracy 
      ind <- which.max(accus)
      text("+", x=ro$specificities[ind], y=ro$sensitivities[ind], pos = 4, cex=1.5, col="red")
      text(paste(c("cutoff = ", as.character(round(ro$thresholds[ind],3)), " (Acc = ", as.character(round(100*accus[ind])/100), ")"), collapse=""), x=ro$specificities[ind], y=ro$sensitivities[ind]-0.05, pos = 4, col="red")    
      cutoff <- ro$thresholds[ind]
    }
    probabilities_cutoff05 <- as.factor(as.numeric(probabilities>=cutoff))
    acc <- sum(probabilities_cutoff05 == Y_test)/length(Y_test)  
    sens <- length(which(Y_test==1 & as.factor(as.numeric(probabilities>=cutoff))==1))/length(Y_test[which(Y_test==1)])
    spec <- length(which(Y_test==0 & as.factor(as.numeric(probabilities>=cutoff))==0))/length(Y_test[which(Y_test==0)])
    
    text(paste("AUC =", round(ro$auc, 3)), x=0.3, y=0.6-2/25, pos = 4, cex=0.6)
    text(paste("Spec at 75% Sens =",round(100*Sens4Spec(ro, 0.75), 1), "%"), x=0.3, y=0.6-3/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at 90% Sens =",round(100*Sens4Spec(ro,  0.9), 1), "%"), x=0.3, y=0.6-4/25, pos = 4, cex=0.6, col="blue")
    text(paste("Best cutoff = ", as.character(round(cutoff,3))), x=0.3, y=0.6-5/25, pos = 4, cex=0.6, col="blue")
    text(paste("Acc at cutoff ", as.character(round(cutoff,3)), " = ",round(100*acc, 1), "%"), x=0.3, y=0.6-6/25, pos = 4, cex=0.6, col="blue")
    text(paste("Sens at cutoff ", as.character(round(cutoff,3)), " = ",round(100*sens, 1), "%"), x=0.3, y=0.6-7/25, pos = 4, cex=0.6, col="blue")
    text(paste("Spec at cutoff ", as.character(round(cutoff,3)), " = ",round(100*spec, 1), "%"), x=0.3, y=0.6-8/25, pos = 4, cex=0.6, col="blue")
    text(subt, x=0.5, y=0, pos=3)
    
    return(cutoff)
  }
}

crossvalidation <- function(predNames = c(), targetname, X, errorfun="aucerror" , seedj = 1, smoo = FALSE, randomtest = FALSE, doScale = FALSE) {
  if (length(predNames)<1) return(c())
  
  set.seed(seedj)
  
  if (randomtest) {
    X[,targetname] <- sample(X[,targetname])
  }
  
  set.seed(seedj) # fuer reproduzierbare Zufallszahlen
  
  modelsAndErrors <- list()
  
  # Bootstrap-Validation bs
  for (i in 1:k) {
    # generiere bs-trainingsset
    bsind <- sample(1:dim(X)[1], dim(X)[1], replace=TRUE)
    if (nlevels(X[,targetname])>0) {
      bsind <- c()
      for (nl in 1:nlevels(X[,targetname])) {
        ind_tmp <- which(X[,targetname]==(levels(X[,targetname])[nl]))
        bsind <- c(bsind, sample(ind_tmp, length(ind_tmp), replace=TRUE))
      }      
    }
    #generiere oob-testset
    oob <- c(1:dim(X)[1])[is.na(match(1:dim(X)[1], bsind))]
    
    trainingsamples <- X[bsind,predNames]
    traininglables <- as.factor(X[bsind,targetname])
    testsamples <- X[oob,predNames]
    testlables <- as.factor(X[oob,targetname])
    
    # feature ranking
    if (n<length(predNames)) {
      importances <- c()
      for (j in 1:100) {
        fitRF <- randomForest(y=as.factor(traininglables), x=trainingsamples,ntree=100, do.trace=F, importance=F, family="binomial")
        if (j==1) {
          importances <- fitRF$importance
        } else {
          importances <- importances + fitRF$importance
        }
      }
      predNames20 <- rownames(importances)[order(importances, decreasing=TRUE)][1:n]
    } else {
      predNames20 <- predNames 
    }
    
    # Exhaustive Search
    for (ind in allInd) {
      predNames_ <- predNames20[ind]
      show(ind)
      markerError <- errorbox_compartment(trainingsamples, traininglables, targetname, cbind(testsamples, testlables), errorfun="prederror", doScale=doScale)
      
      # fuege Modellerror hinzu
      modelsAndErrors <- append(modelsAndErrors, list(id=paste(sort(predNames_), collapse=""), compartments=predNames_, crossValidation=markerError))
    }
    
    
  }
}

compartments2Integercode <- function(compartments, predNames) {
  return(as.real(paste(formatC(which(predNames %in% compartments), width=2, flag="0"), collapse="")))
}

noverk <- function(n,k) {
  return (factorial(n)/(factorial(k)*factorial(n-k)))
}

howmany <- function(n, k) {
  return (sum(sapply(1:k, function(x){noverk(n,x)})))
}

# Leave One Out Cross Correlation: calculates n correlations between vec1 and vec2. Everytime, one item
# of vec1 and vec2 is left out.
loocc <- function(vec1, vec2) {
  cors <- c()
  for (i in 1:length(vec1)) {
    cors <- c(cors, cor(vec1[-i], vec2[-i]))
  }
  return (cors)
}

# weight can be c("unweighted", "equal", "squared")
kappam.light.own <- function (ratings, weight=c("unweighted", "equal", "squared"))
{
  ratings <- as.matrix(na.omit(ratings))
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
    if ((i == 1) & (j == (i + 1)))
      kappas <- kappa2(ratings[, c(i, j)], weight)$value
    else kappas <- c(kappas, kappa2(ratings[, c(i, j)], weight)$value)
  }
  value <- mean(kappas, na.rm=TRUE)
  lev <- levels(as.factor(ratings))
  levlen <- length(levels(as.factor(ratings)))
  if (levlen==1) {
    rval <- structure(list(method = "Light's Kappa for m Raters",
                           subjects = ns, raters = nr, weight = NA, irr.name = "Kappa", value = value,
                           stat.name = "z", statistic = NA, p.value = NA), class = "irrlist")
    return(rval)
  }
  for (nri in 1:(nr - 1)) for (nrj in (nri + 1):nr) {
    for (i in 1:levlen) for (j in 1:levlen) {
      if (i != j) {
        r1i <- sum(ratings[, nri] == lev[i])
        r2j <- sum(ratings[, nrj] == lev[j])
        if (!exists("dis"))
          dis <- r1i * r2j
        else dis <- c(dis, r1i * r2j)
      }
    }
    if (!exists("disrater"))
      disrater <- sum(dis)
    else disrater <- c(disrater, sum(dis))
    rm(dis)
  }
  B <- length(disrater) * prod(disrater)
  chanceP <- 1 - B/ns^(choose(nr, 2) * 2)
  varkappa <- chanceP/(ns * (1 - chanceP))
  SEkappa <- sqrt(varkappa)
  u <- value/SEkappa
  p.value <- 2 * (1 - pnorm(abs(u)))
  rval <- structure(list(method = "Light's Kappa for m Raters",
                         subjects = ns, raters = nr, weight = weight, irr.name = "Kappa", value = value,
                         stat.name = "z", statistic = u, p.value = p.value), class = "irrlist")
  return(rval)
}


AC1m.own <- function (ratings)
{
  ratings <- as.matrix(na.omit(ratings))
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  kappas <- c()
  for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
    lev <- unique(c(levels(factor(ratings[, i])), levels(factor(ratings[, j]))))
    if (length(lev)==1) {
      kappas <- c(kappas, 1)
    } else {
      ratingsi_tmp <- factor(ratings[, i], levels=lev)
      ratingsj_tmp <- factor(ratings[, j], levels=lev)
      kappas <- c(kappas, AC1(table(ratingsi_tmp, ratingsj_tmp),conflev=0.95,N=Inf,print=FALSE)$ac1)
    }
  }
  value <- mean(kappas, na.rm=TRUE)
  lev <- levels(as.factor(ratings))
  levlen <- length(levels(as.factor(ratings)))
  if (levlen==1) {
    rval <- structure(list(method = "AC1 for m Raters",
                           subjects = ns, raters = nr, irr.name = "AC1", value = value,
                           stat.name = "z", statistic = NA, p.value = NA), class = "irrlist")
    return(rval)
  }
  for (nri in 1:(nr - 1)) for (nrj in (nri + 1):nr) {
    for (i in 1:levlen) for (j in 1:levlen) {
      if (i != j) {
        r1i <- sum(ratings[, nri] == lev[i])
        r2j <- sum(ratings[, nrj] == lev[j])
        if (!exists("dis"))
          dis <- r1i * r2j
        else dis <- c(dis, r1i * r2j)
      }
    }
    if (!exists("disrater"))
      disrater <- sum(dis)
    else disrater <- c(disrater, sum(dis))
    rm(dis)
  }
  B <- length(disrater) * prod(disrater)
  chanceP <- 1 - B/ns^(choose(nr, 2) * 2)
  varkappa <- chanceP/(ns * (1 - chanceP))
  SEkappa <- sqrt(varkappa)
  u <- value/SEkappa
  p.value <- 2 * (1 - pnorm(abs(u)))
  rval <- structure(list(method = "AC1 for m Raters",
                         subjects = ns, raters = nr, irr.name = "AC1", value = value,
                         stat.name = "z", statistic = u, p.value = p.value), class = "irrlist")
  return(rval)
}

# AC1 statistic for 2 raters special case
# table = k x k table which represents table(rater1,rater2), must have equal number of rows and columns
# N = population size which will be stick in standard error correction,
# N = Inf is no correction.
# conflev = Confidence Level associated with the confidence interval (0.95 is the default value)
# Gwet KL. Computing inter-rater reliability and its variance in the presence
# of high agreement. Br J Math Stat Psychol. 2008 ;61(Pt 1):29?48.
AC1 <- function(table,conflev=0.95,N=Inf,print=TRUE){
  if(dim(table)[1] != dim(table)[2]){
    stop('The table should have the same number of rows and columns!')
  }
  n <- sum(table)
  f <- n/N
  pa <- sum(diag(table))/n # formula 18
  q <- ncol(table) # number of categories
  pkk <- diag(table)/n
  pak <- sapply(1:q,function(i)sum(table[i,]))/n
  pbk <- sapply(1:q,function(i)sum(table[,i]))/n
  pik <- (pak + pbk)/2
  pegama <- (sum(pik*(1-pik)))/(q-1)
  gama <- (pa - pegama)/(1 - pegama) # AC1 statistics
  # 2 raters special case variance
  pkl <- table/n
  soma <- 0;
  for(k in 1:q){
    for(l in 1:q){
      soma <- soma + (pkl[k,l]*((1-(pik[k]+pik[l])/2)^2))
    }
  }
  vgama <- ((1-f)/(n*(1-pegama)^2)) * (pa*(1-pa) -
                                         4*(1-gama)*((1/(q-1))*sum(pkk*(1-pik)) - pa*pegama) + 4*((1-gama)^2) *
                                         ((1/((q-1)^2))*soma - pegama^2))
  epgama <- sqrt(vgama)# AC1 standard error
  lcb <- max(0,gama - epgama*qnorm(1-(1-conflev)/2,0,1)) # lower confidence bound
  ucb <- min(1,gama + epgama*qnorm(1-(1-conflev)/2,0,1)) # upper confidence bound
  if(print==TRUE){
    cat('Raw agreement:',pa,'Chance-independent agreement:',pegama,'\n')
    cat('Agreement coeficient (AC1):',gama,'AC1 standard
        error:',epgama,'\n')
    cat(conflev*100,'% Confidence Interval (AC1): (',lcb,',',ucb,')\n')
  }
  invisible(c(pa,pegama,gama,epgama,lcb,ucb))
  gwet <- c()
  gwet$raw.agreement <- pa
  gwet$chance.indep.agreement <- pegama
  gwet$ac1 <- gama
  gwet$ac1.stderror <- epgama
  gwet$lowerCI <- lcb
  gwet$upperCI <- ucb
  return (gwet)
}


is.constant <- function(x, na.rm=FALSE) {
  if (na.rm) {
    x <- x[which(!is.na(x))]
  } 
  if (length(x)==0) {
    return (TRUE)
  } else {
    return (all(x==x[1]))
  }
}

is.almost.constant <- function(x, t, na.rm=FALSE) {
  if (na.rm) {
    x <- x[which(!is.na(x))]
  } 
  if (length(x)==0) {
    return (TRUE)
  } else {
    if (any(is.na(x))) {
      return (NA)
    } else {
      s <- summary(as.factor(x))
      return (s[which.max(s)]/sum(s) > 1-t)
    }
  }
}

# calculate the observed Mutual Information between two vectors
mutual_information_from_table <- function(cm, unit="log2") {
  out <- c()
  out$confusion <- cm
  out$entropy_v1 <- entropy.plugin(rowSums(cm), unit=unit)
  out$entropy_v2 <- entropy.plugin(colSums(cm), unit=unit)
  out$mi <- mi.plugin(cm, unit=unit)
  out$mi_normalized <- out$mi / max(out$entropy_v1, out$entropy_v2)  
  
  return (out)
}

# calculate the observed Mutual Information between two vectors
mutual_information <- function(v1, v2, unit) {
  #if (is.numeric(v1)) v1 = as.factor(v1>median(v1, na.rm=TRUE))
  #if (is.numeric(v2)) v2 = as.factor(v2>median(v2, na.rm=TRUE))
  cm <- table(v1, v2)
  cm <- cm/sum(cm)
  return (mutual_information_from_table(cm, unit))
}


MyScale <- function(X, ...) {
  out <- scale(X, ...)
  out <- sapply(1:dim(X)[2], function(x) {(if (all(is.na(out[,x]))) rep(0, length(out[,x])) else out[,x])})
  return(out)
}


entropy <- function(X, epsilon=0) {
  X <- (X+epsilon)/(sum(X)+epsilon)
  H <- -sum(sapply(X, function(x) {(if (!is.na(x) && x>0) x*log2(x) else 0)}))
  return(H)
}

cond_entropy <- function(ctable, epsilon=0) {
  cond_P <- (ctable+epsilon)/(t(matrix(rep(colSums(ctable),dim(ctable)[1]), nrow=dim(ctable)[2]))+epsilon)
  P <- colSums(ctable/sum(ctable));
  H <- 0;
  for (i in 1:dim(ctable)[2]) {
    H <- H + P[i]*entropy(cond_P[,i])
  }
  return(as.numeric(H))
}


mutual_information2 <- function(v1, v2, epsilon=0) {
  #if (is.numeric(v1)) v1 = as.factor(v1>median(v1, na.rm=TRUE))
  #if (is.numeric(v2)) v2 = as.factor(v2>median(v2, na.rm=TRUE))
  
  ctable = table(v1, v2)
  out <- c()
  out$confusion <- ctable
  out$entropy_v1 <- entropy(rowSums(ctable), epsilon)
  out$entropy_v2 <- entropy(colSums(ctable), epsilon)
  out$mi <- entropy(rowSums(ctable), epsilon) - cond_entropy(table(as.numeric(v1), as.numeric(v2)), epsilon)
  out$mi_normalized <- out$mi / max(out$entropy_v1, out$entropy_v2)
  
  return(out)
}

my_pairs <- function(X, col=null, ...) {
  p.old <- par(no.readonly=TRUE)
  par(mfrow=c(dim(X)[2],dim(X)[2]), mar=c(0.5,0.5,0.5,0.5), oma=c(3,3,3,3))
  for (i in 1:dim(X)[2]) {
    for (j in 1:dim(X)[2]) {
      if (i==j) {
        plot(c(0,1), c(0,1), type="n", xlab="", ylab="", xaxt="n", yaxt="n", ...)
        text(0.5, 0.5, colnames(X)[i], cex=0.9)
      } else {
        if (length(unique(X[,i]))>4 || length(unique(X[,j]))>4) {
          plot(X[,i], X[,j], col=col, xlab="", ylab="", xaxt="n", yaxt="n", ...)
        } else {
          f <- as.data.frame(table(X[,i],X[,j]), responseName = "Freq", stringsAsFactors = TRUE)
          col2 <- c()
          xs <- c()
          ys <- c()
          cexs <- c()
          for (k in 1:dim(f)[1]) {
            if (f$Freq[k] > 0) {
              xs <- c(xs, as.numeric(as.character(f[k,1])))
              ys <- c(ys, as.numeric(as.character(f[k,2])))
              cexs <- c(cexs, as.numeric(as.character(f[k,3])))
              cols <- as.data.frame(table(c=col[which(X[,i]==f[k,1] & X[,j]==f[k,2])]))
              col2 <- c(col2, as.numeric(as.character(cols$c[which(cols$Freq==max(cols$Freq))[1]])))
            }
          }
          cexs <- 10*cexs/max(cexs)
          plot(xs, ys, col=col2, xlab="", ylab="", xaxt="n", yaxt="n", cex=cexs, ...)
        }
        if (i==1 && j%%2==0) {
          axis(3)
        } else if (i==dim(X)[2] && j%%2==1) {
          axis(1)
        }
        if (j==1 && i%%2==0) {
          axis(2)
        } else if (j==dim(X)[2] && i%%2==1) {
          axis(4)
        }
      }
    }
  }
  par(p.old)
}

findModelIndex <- function(predNames, predNames20, bsindAll, exact=TRUE) {
  for (ind_tmp in 1:length(bsindAll)) {
    preds <- predNames20[unlist(bsindAll[ind_tmp])]
    if ((!exact || length(preds)==length(predNames)) && all(predNames %in% preds)) {
      return(ind_tmp)
    }
  }
  return (NA)
}

# Lopez-Paz, Hennig, Sch?lkopf: The Randomized Dependence Coefficient. Stat ML 2013.
rdc <- function(x, y, k, s) {
  if (is.na(var(x, na.rm=TRUE)) || var(x, na.rm=TRUE)==0 || is.na(var(x, na.rm=TRUE)) || var(y, na.rm=TRUE)==0) {
    return (NA)
  }
  x <- cbind(apply(as.matrix(x), 2, function(u) ecdf(u)(u)), 1)
  y <- cbind(apply(as.matrix(y), 2, function(u) ecdf(u)(u)), 1)
  wx <- matrix(rnorm(ncol(x)*k, 0, s), ncol(x), k)
  wy <- matrix(rnorm(ncol(y)*k, 0, s), ncol(y), k)
  cancor(cbind(cos(x%*%wx), sin(x%*%wx)), cbind(cos(y%*%wy), sin(y%*%wy)))$cor[1]
}

cor.mtest <- function (mat, conf.level = 0.95, method="spearman") 
{
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level, method=method)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      if (method=="pearson") {
        lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
        uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
      }
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

plotCorrelations <- function (data, maintitle = "", method="spearman", ...) {
  require(corrplot)
  cor_jnk=cor(data, use="pairwise.complete.obs", method = method)
  cor_jnk[which(is.na(cor_jnk))] <- 0
  tmp <- cor.mtest(data, method=method)
  corrplot(cor_jnk, method="circle", tl.pos="lt", type="upper",
           p.mat = tmp[[1]], pch.col="red", sig.level = 0.05,
           tl.col="black", tl.srt=45, addCoef.col=c("black", "white")[as.numeric(abs(cor_jnk[upper.tri(cor_jnk, diag = TRUE)])>0.6)+1], addCoefasPercent = TRUE, main="", ...)
  title(main = maintitle, line=3, cex.main = 2, ...)
  title(sub = paste("Spearman Rank Correlation r, n =", nrow(data), collapse=""), line = 0, cex.sub=2, ...)
  
  result <- c()
  result$cor <- cor_jnk
  result$cortest <- tmp
  return(result)
}

