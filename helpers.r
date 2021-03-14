###
### Peter J. Schueffler
###
### Date: 2015-10-29
###

require(doParallel)
require(sampling)
require(pROC)

##
## LSF HPC management
##

jobsnmax = function () 
{
  defnmax <- 288
  cmd <- "busers"
  out <- system(cmd, intern = TRUE, ignore.stderr = FALSE, 
                wait = TRUE, input = NULL)
  if (grepl("Cannot connect to LSF. Please wait ...", out[1])) {
    return(defnmax)
  }
  strsplit(out[1], " +")[[1]]
  out <- strsplit(out[2], " +")[[1]]
  return(tryCatch(as.numeric(out[3]), warning = function(e) {
    defnmax
  }, error = function(e) {
    defnmax
  }))
}

jobsFinished = function (skip = 0) 
{
  cmd <- "busers"
  out <- system(cmd, intern = TRUE, ignore.stderr = FALSE, 
                wait = TRUE, input = NULL)
  if (grepl("Cannot connect to LSF. Please wait ...", out[1])) {
    return(FALSE)
  }
  strsplit(out[1], " +")[[1]]
  out <- strsplit(out[2], " +")[[1]]
  return(out[5] == "0" && !is.na(as.numeric(out[6])) && as.numeric(out[6]) <= 
           skip)
}

waitUntilJobsFinished = function (skip = 0) 
{
  Sys.sleep(5)
  while (!jobsFinished(skip)) {
    Sys.sleep(60)
  }
  return(TRUE)
}

delete_successful_lsf = function (foldername = ".") 
{
  lsf_files <- dir(foldername, pattern = "lsf.*")
  if (length(lsf_files) > 0) {
    for (i in 1:length(lsf_files)) {
      words = read.table(lsf_files[i], sep = "\n", stringsAsFactors = FALSE, 
                         nrows = 30)
      if (any(words == "Successfully completed.")) {
        unlink(lsf_files[i])
      }
    }
  }
}

##
## To create all 1-5 model combinations
##

CreateIndices_fast = function (n, p) 
{
  cmd <- paste(c("sapply(1:(", as.character(n), "-", as.character(p - 
                                                                    1), "), function(v", as.character(1), "){PLACEHOLDER})"), 
               collapse = "")
  if (p > 1) {
    for (i in 2:p) {
      cmd <- sub("PLACEHOLDER", paste(c("sapply((v", as.character(i - 
                                                                    1), "+1):(", as.character(n), "-", as.character(p - 
                                                                                                                      i), "), function(v", as.character(i), "){PLACEHOLDER})"), 
                                      collapse = ""), cmd)
    }
  }
  cmd <- sub("PLACEHOLDER", paste(c("list(c(", paste(sapply(1:p, 
                                                            function(x) {
                                                              paste(c("v", as.character(x)), collapse = "")
                                                            }), collapse = ","), "))"), collapse = ""), cmd)
  cmd <- paste(c("unlist(", cmd, ")"), collapse = "")
  tmp <- eval(parse(text = cmd))
  out <- lapply(1:(length(tmp)/p), function(x) {
    c(tmp[(p * (x - 1) + 1):(p * x)])
  })
  return(out)
}

CreateIndices = function (n, p) 
{
  out <- list()
  for (i in 1:p) {
    out <- append(out, CreateIndices_fast(n, i))
  }
  return(out)
}

findModelIndex = function (predNames, predNames20, bsindAll, exact = TRUE) 
{
  for (ind_tmp in 1:length(bsindAll)) {
    preds <- predNames20[unlist(bsindAll[ind_tmp])]
    if ((!exact || length(preds) == length(predNames)) && 
        all(predNames %in% preds)) {
      return(ind_tmp)
    }
  }
  return(NA)
}



##
## For Plots
##

shortenNames = function (namelist, reallyshort = FALSE) 
{
  if (reallyshort) {
    return(sapply(namelist, function(x) {
      strsplit(x, "_")[[1]][1]
    }, USE.NAMES = FALSE))
  }
  else {
    return(sapply(namelist, function(x) {
      l = strsplit(x, "_")[[1]]
      if (length(l) > 1) {
        paste(c(l[1], substr(l[2], 1, 1)), collapse = "_")
      } else {
        l[1]
      }
    }, USE.NAMES = FALSE))
  }
}
  

predHist2 = function (allpreds, tit = "", short.names = FALSE, cex.names = 0.6, 
                      mai = c(1.02, 1.35, 0.82, 0.42), n = -1, normalized = FALSE, 
                      doPlot = TRUE, barcol = "orange", ...) 
{
  his <- unique(unlist(allpreds))
  his <- sapply(his, function(x) {
    sum(unlist(allpreds) == x)
  })
  his <- sort(his)
  if (n > 0) {
    his <- his/n
  }
  if (normalized) {
    his <- his/sum(his)
  }
  if (short.names) {
    names(his) <- shortenNames(names(his), reallyshort = TRUE)
  }
  if (doPlot) {
    par(mai = mai)
    if (n > 0) {
      his <- 100 * his
    }
    barplot(his, horiz = TRUE, col = barcol, las = 1, cex.names = cex.names, 
            main = tit, ...)
    par(mai = c(1.02, 0.82, 0.82, 0.42))
  }
  return(his)
}

errormeasureAUC = function (y, y_hat, doPlot = FALSE, smooth = FALSE, getROC = FALSE, ...) 
{
  roc_ = roc(y, y_hat, smooth = smooth)
  if (doPlot) {
    plot(roc_, mar = c(5.1, 4.1, 4.2, 2.1), ...)
  }
  if (getROC) {
    return(roc_)
  }
  else {
    return(auc(roc(y, y_hat)))
  }
}

myboxplot = function (x, ci = 0.95, ...) 
{
  tryCatch({
    myboxplot1(x, ci, ...)
  }, error = function(y) {
    myboxplot2(x, ci, ...)
  })
}

myboxplot1 = function (x, ci = 0.95, ...) 
{
  out <- boxplot(x, ...)
  medians <- apply(x, 2, function(y) median(y, na.rm = TRUE))
  cilower <- apply(x, 2, function(y) {
    mycimedian(y)[2]
  })
  ciupper <- apply(x, 2, function(y) {
    mycimedian(y)[3]
  })
  segments(1:length(medians), cilower, 1:length(medians), ciupper, 
           lwd = 10, col = "lightgrey")
  points(medians, col = "black", cex = 2.5, pch = "-")
  return(out)
}

myboxplot2 = function (x, ci = 0.95, ...) 
{
  out <- boxplot(x, ...)
  if (typeof(x) != "list") {
    x <- list(x)
  }
  medians <- sapply(x, function(y) median(y, na.rm = TRUE))
  cilower <- sapply(x, function(y) {
    mycimedian(y)[2]
  })
  ciupper <- sapply(x, function(y) {
    mycimedian(y)[3]
  })
  segments(1:length(medians), cilower, 1:length(medians), ciupper, 
           lwd = 10, col = "lightgrey")
  points(medians, col = "black", cex = 2.5, pch = "-")
  return(out)
}

mycimedian = function (x) 
{
  n = length(x)
  x_sorted = sort(x)
  out = c(median(x), x_sorted[floor((n + 1)/2 - 1.96 * sqrt(n)/2)], 
          x_sorted[ceiling((n + 1)/2 + 1.96 * sqrt(n)/2)])
  return(out)
}

avg.roc = function (rocs) 
{
  rocs <- rocs[which(!is.na(rocs))]
  ro <- rocs[[1]]
  xs <- (0:100)/100
  ys <- xs - xs
  for (i in 1:length(xs)) {
    ys[i] <- median(sapply(1:length(rocs), function(x) {
      get.sens.at(xs[i], rocs[[x]])
    }))
  }
  ys[1] = 0
  ys[length(xs)] = 1
  ro$specificities <- xs
  ro$sensitivities <- ys
  ro$auc = 0
  for (i in 2:length(xs)) {
    ro$auc = ro$auc + (1/length(xs) * (ys[i] + ys[i - 1])/2)
  }
  return(ro)
}

get.sens.at = function (x, ro) 
{
  x1 <- ro$specificities
  y1 <- ro$sensitivities
  m <- 0
  if (is.element(x, x1)) {
    m <- max(y1[which(x1 == x)])
  }
  else {
    ya <- x1[rev(which(x1 < x))[1]]
    yb <- x1[(which(x1 > x))[1]]
    a <- y1[rev(which(x1 < x))[1]]
    b <- y1[(which(x1 > x))[1]]
    m <- a + (b - a) * (x - ya)/(yb - ya)
  }
  return(m)
}

##
## General
##

removeNas = function (data, predNames) 
{
  nas <- which(is.na(data[, c(predNames)]), arr.ind = TRUE)
  if (length(nas) > 0) {
    if (length(dim(nas)) == 0) {
      data <- data[-unique(nas), ]
    }
    else if (length(nas[, 1]) > 0) {
      data <- data[-unique(nas[, 1]), ]
    }
  }
  return(data)
}