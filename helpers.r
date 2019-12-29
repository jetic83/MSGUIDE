################################################################################
##
##  Helper Functions
##
##  Thomas J. Fuchs
##  ETH Zurich, Institute for Computational Science
##  thomas.fuchs@inf.ethz.ch
##
################################################################################


s <- summary
d.f <- data.frame
g.o <- graphics.off



mb <- function(x)
{
   return(object.size(x)/1024^2)
}


stamp <- function()
{
   z <- format(Sys.time(), "%c")
   z <- gsub(":", "-", z)
   z <- gsub(" ", "_", z)
   return(z)
}


################################################################
### List Handling
################################################################


list.colapse <- function(datalist) {
   newList <- list();
   for (i in 1:length(datalist)) {
      if((length(newList)==0)||(newList[length(newList)] != datalist[i])) {
         newList <- c(newList, datalist[i]);
      }
   }
   # return value:
   newList
}


list.normalize <- function(l) {
	if (length(which(l<0))>0) {
		print("ERROR in function list.normalize \n NEGATIVE LIST ENTRIES!");
	} else {
		s <- sum(l);
		l <- l/s;
	}
	# return value:
	l
}

list.negToZero <- function(l) {
	w <- which(l < 0);
	for (i in 1:length(w)) {
		l[w[i]] <- 0;
	}
	# return value:
	l
}



################################################################
### Data Frame Handling
################################################################

numeric.data.frame <- function(dFrame)
{
   numFrame <- dFrame
   for (i in 1:ncol(dFrame))
   {
      if (is.factor(dFrame[,i]))
      {
         numFrame[,i] <- as.numeric(dFrame[,i]) - 1
      } else {
         if (!is.numeric(dFrame[,i]))
         {
            # Sanity Check
            stop(paste("Column",i,"is neither numeric nor a factor! It is a ", class(dFrame[,i])))
         }
      }
   }
   rownames(numFrame) <- rownames(dFrame)
   names(numFrame) <- names(dFrame)

   return(numFrame)
}



factor.data.frame <- function(dFrame, fLevels)
{
   nFrame <- dFrame
   for (i in 1:ncol(nFrame))
   {
      nFrame[,i] <- as.factor(nFrame[,i])
      levels(nFrame[,i]) <- fLevels
   }
   return(nFrame)
}

character.data.frame <- function(dFrame)
{
   nFrame <- dFrame
   for (i in 1:ncol(nFrame))
   {
      nFrame[,i] <- as.character(nFrame[,i])
   }
   return(nFrame)
}



sort.data.frame <- function(x, key, ...) {
    if (missing(key)) {
        rn <- rownames(x)
        if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
        x[order(rn, ...), , drop=FALSE]
    } else {
        x[do.call("order", c(x[key], ...)), , drop=FALSE]
    }
}

sample.data.frame <- function(df, ...) {
   ids <- sample(c(1:nrow(df)), ...);
   select.rows(df, ids)
}



types.data.frame <- function(X)
{
   types <- matrix(F,nrow=ncol(X),ncol=3)
   colnames(types) <- c("Factor", "Numeric", "Character")

   for (i in 1:ncol(X))
   {
      types[i,1] <- is.factor(X[,i])
      types[i,2] <- is.numeric(X[,i])
      types[i,3] <- is.character(X[,i])
      if (sum(types[i,]) != 1)
      {
         stop(paste("Unknown type in column",i))
      }
   }
   return(data.frame(types))
}


add.rows.to.data.frame <- function(frame1,frame2) {
   indToAdd <- c((nrow(frame1)+1):(nrow(frame1)+nrow(frame2)))
   frame1[indToAdd,] <- frame2
   return(frame1)
}


data.frame.standardize <- function(X)
{
   for (i in 1:ncol(X))
   {
      if (is.numeric(X[,i]))
      {
         if (max(X[,i],na.rm=TRUE) != 0)
         {
            X[,i] <- X[,i] / max(X[,i],na.rm=TRUE)
         }
      }
   }
   return(X)
}


################################################################
### Plotting
################################################################

# http://tolstoy.newcastle.edu.au/R/help/06/04/26017.html
image.matrix <-
    function (x, rlab = if (is.null(rownames(x))) 1:nrow(x) else rownames(x),
              clab = if (is.null(colnames(x))) 1:ncol(x) else colnames(x),
              cex.row=0.7, cex.col=0.7, main = deparse(substitute(x)),
              numbers=TRUE, ...)
{
        op <- par(mgp=c(2, .3, 0))
        #on.exit(par(op))
        nr <- nrow(x)
        nc <- ncol(x)
        image(1:nc, 1:nr, t(x)[, nr:1], axes = FALSE, xlab = "",
          ylab = "", main = main, ...)
        axis(2, 1:nr, rev(rlab), cex.axis=cex.col, tick=FALSE, las=2)
        axis(1, 1:nc, clab, cex.axis=cex.row, tick=FALSE, las=2)

        if(numbers) {
          #cmt <- t(cm)
          for (r in 1:nrow(cm)) {
            for (c in 1:ncol(cm)) {
              #text(x=c, y=r, format(cm[r,c], digits=1, nsmall=2))
              #text(x=ncol(cm)-c+1, y=r, format(cm[r,c], digits=1, nsmall=2))
              text(x=c, y=nrow(cm)-r+1, format(cm[r,c], digits=1, nsmall=2))
            }
          }
        }
        invisible(x)
}


smoothPairs <- function(X, ...)
{
   library(geneplotter)
   par(mfrow=c(1,1))
   pairs(X, panel=function(...) {par(new=TRUE);smoothScatter(..., nrpoints=0)})

   par(oldpar)
}





################################################################
### NA Handling
################################################################

percentageMissing <- function(X)
{
   return(sum(is.na(X)) / (dim(X)[1]*dim(X)[2]))
}


replaceNAsMean <- function(X)
{
   means <- apply(X, 2, mean, na.rm=TRUE)
   for (i in 1:ncol(X))
   {
      if (sum(is.na(X[,i]))>0)
      {
         X[which(is.na(X[,i])),i] <- means[i]
      }
   }
   return(X)
}



replaceNAsMedian <- function(X)
{
   medians <- apply(X, 2, median, na.rm=TRUE)
   for (i in 1:ncol(X))
   {
      if (sum(is.na(X[,i]))>0)
      {
         X[which(is.na(X[,i])),i] <- medians[i]
      }
   }
   return(X)
}


df.removeNaRows <- function(X)
{
   ind <- vector()
   for (i in 1:ncol(X))
   {
      ind <- c(ind,which(is.na(X[,i])))
   }
   ind <- unique(ind)
   return(X[-ind,])
}


removeNaRows <- function(X1,X2)
{
   aind <- which(is.na(X1), arr.ind=TRUE)
   nullRows <- sort(unique(aind[,1]))
   return(list(X1=X1[-nullRows,],
               X2=X2[-nullRows,]))
}


data.frame.fillUpNas <- function(X)
{
   types <- types.data.frame(X)

   # Replace numircal NAs with mean:
   indNum <- which(types$Numeric==TRUE)
   if (length(indNum)>0)
   {
      #X[,indNum] <- replaceNAsMean(X[,indNum])
      X[,indNum] <- replaceNAsMedian(X[,indNum])
   }

   # Replace factor NAs with regression from the numeric ones:
   indFactor <- which(types$Factor==TRUE)
   if (length(indFactor)>0)
   {
      for (i in 1:length(indFactor))
      {
         if (sum(is.na(X[,indFactor[i]])) > 0)
         {
            X[,indFactor[i]] <- fillMissingBinaryFactor(X[,indFactor[i]],X[,indNum])
         }
      }
   }
   
   return(X)
}



fillMissingBinaryFactor <- function(f1,Xnum)
{
   indNA <- which(is.na(f1))
   trainX <- Xnum[-indNA,]
   testX <- Xnum[indNA,]
   fit <- glm(f1[-indNA]~., data=trainX, family=binomial())
   f1_hat <- predict(fit,testX)
   f1_hat.bin <- (f1_hat>0)
   f1[indNA][f1_hat.bin] <- levels(f1)[1]
   f1[indNA][!f1_hat.bin] <- levels(f1)[2]

   return(f1)
}


################################################################
### Cross Validation
################################################################

getKFoldBlocks <- function(sampleSize, folds)
{
   indices <- 1:sampleSize
   indices <- sample(indices)
   blockSize <- as.integer(sampleSize/folds)
   blocks <- rep(blockSize, folds)
   
   #distribute the rest of the instances to the first blocks
   mod <- sampleSize%%folds
   blocks[1:mod] <- blocks[1:mod] + 1
   
   foldList <- list()
   for (i in 1:folds)
   {
      foldList[[i]] <- indices[c( (sum(blocks[1:i-1])+1) : sum(blocks[1:i]) )]
   }
   
   return(foldList)
}



################################################################
### Miscelanous
################################################################



dens2d <- function(x, nx = 20, ny = 20, margin = 0.05, h = 1)
{
   xrange <- max(x[, 1]) - min(x[, 1])
   yrange <- max(x[, 2]) - min(x[, 2])
   xmin <- min(x[, 1]) - xrange * margin
   xmax <- max(x[, 1]) + xrange * margin
   ymin <- min(x[, 2]) - yrange * margin
   ymax <- max(x[, 2]) + yrange * margin

   xstep <- (xmax - xmin)/(nx - 1)
   ystep <- (ymax - ymin)/(ny - 1)
   xx <- xmin + (0:(nx - 1)) * xstep
   yy <- ymin + (0:(ny - 1)) * ystep
   g <- matrix(0, ncol = nx, nrow = ny)
   n <- dim(x)[[1]]
   for(i in 1:n) {
      coefx <- dnorm(xx - x[i, 1], mean = 0, sd = h)
      coefy <- dnorm(yy - x[i, 2], mean = 0, sd = h)
      g <- g + coefx %*% t(coefy)/n
   }
   return(list(x = xx, y = yy, z = g))
}



table.ordered <- function(dataTable, orderList) {
   frequencies <- table(dataTable);
   newTable <- list();
   for(caption in orderList) {
      newTable <- cbind(newTable, frequencies[caption]);
   }
   names(newTable) <- orderList;
   # return value:
   newTable
}

select.columns <- function(dataTable, selection) {
   newTable <- dataTable[selection[1]];
   for(i in 2:length(selection)) {
      newTable <- cbind(newTable, dataTable[selection[i]]);
   }
   # return value:
   newTable
}

select.rows <- function(dataTable, selection){
   newTable <- dataTable[selection[1],];
   for(i in 2:length(selection)) {
      newTable <- rbind(newTable, dataTable[selection[i],]);
   }
   # return value:
   newTable
}


norm <- function(x){ 
   sqrt(x%*%x)
}

vectors.angle <- function(a,b){
   acos(((a%*%b)/(norm(a)*norm(b)))-1e-7)
}

radToDegree <- function(rad) {
   (180/pi)*rad   
}

r2d <- function(rad) {
   radToDegree(rad)
}

con <- function(string1, string2) {
   paste(string1, string2, sep="")
}



# freq displays a simple frequency table for a vector, or steps through
# the columns of a data frame or matrix and returns the result of tabulate()
freq <- function(x,labels,nbins,var.label,maxbins=10,width=8,show.pc=F) {
    if(!missing(x)) {
     if(is.null(dim(x))) {
      # see if there are any NAs
      nna<-sum(is.na(x))
      if(is.numeric(x)) {
       if(missing(labels)) {
        # tabulate drops categories for integers having no observations
        xrange<-range(na.omit(x))
        categories<-seq(xrange[1],xrange[2])
       }
       else categories<-labels
       # another serious kludge - tabulate() always starts at 1
       if(xrange[1] > 1) x<-x-(xrange[1]-1)
      }
      if(is.factor(x)) categories<-levels(x)
      # serious kludge to get around tabulate() shifting 0 values in factors
      if(categories[1] == "0") x<-x+1
      # get the variable label here or it might be clobbered
      if(missing(var.label)) var.label<-deparse(substitute(x))
      if(missing(nbins)) nbins<-length(categories)
      # this effectively clobbers the maxbins limit - try to fix
      if(nbins > maxbins) maxbins<-nbins
      if(length(categories) <= maxbins) {
       # if NAs present, tack on a label
       if(nna) categories<-c(categories,"NA")
       categories<-formatC(ifelse(categories=="","Missing",categories),width=width)
       cat("Frequencies for",var.label,"\n")
       cat(" ",categories,"\n")
       # tabulate barfs with NAs
       freqs<-tabulate(na.omit(x),nbins)
       # tack on the NA count
       if(nna) freqs<-c(freqs,nna)
       cat(" ",formatC(as.character(freqs),width=width),"\n")
       if(show.pc) {
        cat("%",formatC(as.character(round(100*freqs/sum(freqs),1)),width=width),"\n")
       }
       cat("\n")
       names(freqs)<-categories
       invisible(freqs)
      }
      else cat(length(categories),"categories exceeds maximum bins!\n")
     }
     else {
      nfreq<-dim(x)[2]
      freq.list<-rep(list(0),nfreq)
      var.labels<-names(x)
      for(i in 1:nfreq)
       freq.list[[i]]<-freq(x[,i],labels,nbins,var.labels[i],maxbins,width,show.pc)
      names(freq.list)<-names(x)
      invisible(freq.list)
     }
    }
    else cat("Usage: freq(x,labels=NULL,nbins,maxbins=10,width=10,show.pc=F)\n")
}


# Convert all comlums in a data.frame to factors 
# with a given set of levels.
toFactor <- function(dataFrame, levelsList) {
   for (i in 1:ncol(dataFrame))
   {
      dataFrame[[i]] <- factor(dataFrame[[i]], levels=levelsList)
   }
   invisible(dataFrame)
}


#
# Creates two way contingency tables
# from: http://tolstoy.newcastle.edu.au/R/help/01a/1483.html
#
twoWay <- function( x=NA, y=NA, userDefined=NA ){ 

  if (is.na(userDefined)){ 
    result <- chisq.test(table(x,y)) 
  } 
  else{ 
    result <- chisq.test(userDefined) 
  } 

  print (result) 
  observed <-result$observed 
  expected <- result$expected 
  chi.table <- ((observed - expected)^2)/expected 
  row.sum <- apply(observed,1,sum) 
  col.sum <- apply(observed,2,sum) 
  N <- sum(observed) 

  ## put in the marginals and names ... create fullArray 
  fullArray <- cbind(observed,row.sum) 
  fullArray <- rbind(fullArray,c(col.sum,N)) 
  rownames(fullArray) <- c(rownames(observed), "Total") 
  colnames(fullArray) <- c(colnames(observed), "Total") 

  ## make the tables of proportions                       
  proportion <- fullArray/N 
  row.proportion <- fullArray/c(row.sum,N) 
  col.proportion <- t(t(fullArray)/c(col.sum, N)) 

  return(list(fA=fullArray, e=expected, ct=chi.table, p=proportion, 
         rp=row.proportion, cp=col.proportion)) 
}


cout <- function(string) {
   cat(string,"\n")
   if (.Platform$OS.type == "windows") flush.console()
}


writeForWeka <- function(fileName, dataFrame, ...) {
   write.table(dataFrame, fileName, row.names=F, col.names=TRUE, sep=",", quote=TRUE, ...);
}

write.weka <- writeForWeka


replicateClasses <- function(X, Y) {
   countClass1 <- length(which(Y==1))   
   countClass2 <- length(which(Y==2))   
   
   if (countClass1 < countClass2) {
      indToReplicate <- which(Y==1)   
      timesReplicate <- floor(countClass2/countClass1)-1
   } else {
      indToReplicate <- which(Y==2)      
      timesReplicate <- floor(countClass1/countClass2)-1
   }   

   if (timesReplicate > 0) {
      repData <- data.frame(X,Y)
      toRepDat <- repData[indToReplicate,]
      for (i in 1:timesReplicate) {
         repData <- add.rows.to.data.frame(repData,toRepDat)
      }
      
      Y <- repData$Y
      X <- repData[, -c(which(names(repData)=="Y"))]      
   }
   return(list(X=X,Y=Y))
}



BER <- function(truth, prediction)
{
   tt <- table(truth, prediction)
   
   if (max(dim(tt)) > 2)
   {
      stop("The function BER is only defined for two level factors!")
   }
   
   if (dim(tt)[1] != dim(tt)[2])
   {
      BER <- 0.5
   } else {
      a <- tt[1,1]
      b <- tt[1,2]
      c <- tt[2,1]
      d <- tt[2,2]
      BER <- 0.5*(b/(a+b) + c/(c+d))
   }
   return(BER)
}



removeIdenticalCols <- function(X)
{
   numX <- numeric.data.frame(X)
   nonNA  <- apply(!is.na(numX), 2, sum)
   colSum <- apply(numX, 2, sum, na.rm=TRUE)
   indIdentical <- which((colSum%%nonNA)==0)
   return(X[,-indIdentical])
}



simpleInteractions <- function(X)
{
   X2 <- X
   newNames <- dimnames(X)[[2]]
   for (i in 1: ncol(X))
   {
      X2 <- cbind(X2, X[,i:ncol(X)]*X[,i])
      if (i==ncol(X))
      {
         newNames <- c(newNames, paste(dimnames(X)[[2]][i], dimnames(X)[[2]][i], sep=":"))
      } else {
         newNames <- c(newNames, paste(dimnames(X[,i:ncol(X)])[[2]], dimnames(X)[[2]][i], sep=":"))
      }
   }
   dimnames(X2)[[2]] <- newNames
   return(X2)
}


loadExcelTable <- function(fileName, tableName)
{
   require(RODBC)
   channel <- odbcConnectExcel2007(fileName)
   #tables <- sqlTables(channel)
   #tableNames <- substr(tables$TABLE_NAME, 2, nchar(tables$TABLE_NAME)-2)
   O <- sqlFetch(channel, tableName)   
   #drug$name <- sqlColumns(channel, tableNames[i])$COLUMN_NAME[1]
   
   odbcCloseAll()                                                 
   return(O)
}


################################################################
### ROCs
################################################################


getSens <- function(tb)
{
   TP <- tb[1,1]
   TN <- tb[2,2]
   FP <- tb[1,2]
   FN <- tb[2,1]
   return(TP/(TP+FN))
}

getSpec <- function(tb)
{
   TP <- tb[1,1]
   TN <- tb[2,2]
   FP <- tb[1,2]
   FN <- tb[2,1]
   return(TN/(TN+FP))
}



################################################################
### String Manipulation
################################################################

string.trim <- function(string)
{
   string <- gsub('^[[:space:]]+', '', string)
   string <- gsub('[[:space:]]+$', '', string)
   return(string)
}

twoDigits  <- function(txt)
{
   if (nchar(txt) < 2)
   {
      return(paste("0",txt,sep=""))   
   } else {
      return(txt)
   }        
}




################################################################
### Survival
################################################################



reversCensoring <- function(cc)
{
   cc[cc==0] <- 2
   cc[cc==1] <- 0
   cc[cc==2] <- 1
   return(cc)
}



################################################################
### Random Forest
################################################################

rf.importance<- function(forest)
{
   imp <- as.vector(forest$importance)
   names(imp) <- rownames(forest$importance)
   impSorted <- sort(imp, decreas=TRUE)
   X11()
   plot(impSorted)
   names(impSorted[1:50])
   return(impSorted)
}

rf.imp.barplot <- function(forest, ...)
{
   imp <- as.vector(forest$importance)
   names(imp) <- rownames(forest$importance)
   impSorted <- sort(imp, decreas=F)
   X11()
   par(mar=c(5,8,2,2))
   barplot(impSorted, horiz=T, las=2, cex.names=0.7,  col="orange", ...)
}



################################################################
### short cuts
################################################################

df.numeric <- numeric.data.frame
df.factor <- factor.data.frame
df.character <- character.data.frame
df.sort <- sort.data.frame
df.sample <- sample.data.frame
df.classes <- types.data.frame
df.add <- add.rows.to.data.frame
df.standardize <- data.frame.standardize
df.fillUpNas <- data.frame.fillUpNas
df.types <- types.data.frame

