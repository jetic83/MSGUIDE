#bsub -W 08:00 R --file=run_automated.r --no-save --args 1 0 4
#R < run_automated.r --no-save

source("helpers.r")

for (i in 1:length(commandArgs())) {
  show(commandArgs()[i])
}
seedj_tmp <- 1
if (length(commandArgs())>4) {
  seedj_tmp <- as.numeric(commandArgs()[5])
}
isMetaRun <- 0
if (length(commandArgs())>5) {
  isMetaRun <- as.numeric(commandArgs()[6])
}

batchj <- 4
if (length(commandArgs())>6) {
  batchj <- as.numeric(commandArgs()[7])
}

show(seedj_tmp)
show(isMetaRun)

show(batchj)

# Directories
wd <- getwd()
td <- sub("/peschuef/", "/peschuef/work/", wd)
if (!file.exists(td)) dir.create(td, recursive = TRUE)

# Script 1
cmd <- paste(c("bsub R --file=prostate_Script1.r --save --args ", as.character(seedj_tmp), " ", as.character(batchj)), collapse="")
system(cmd, intern = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
waitUntilJobsFinished(1+isMetaRun)
load(".RData")

## Script 2
z <- 1
numberModels <- length(allInd)
bnmax <- jobsnmax()-1-isMetaRun
stepsize <- 50
i <- 1
while(stepsize>1*120) {
  stepsize <- ceiling(numberModels/(i*bnmax))
  i <- i+1
}
save.image()

unlink(paste(c(td, "/*.csv"), collapse=""))
for (j in 1:z) {
  for (i in 1:ceiling(numberModels/stepsize)) {
    tablename <- paste(c(td, "/table", j, "-", i, ".csv"), collapse="")
    cmd <- paste(c("bsub -W 01:00 R --file=prostate_Script2.r --no-save --args ", tablename, " ", j, " ", stepsize*(i-1)+1, " ", min(stepsize*i, numberModels)),  collapse="")
    system(cmd, intern = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
  }
}


waitUntilJobsFinished(1+isMetaRun)


## Script 3
# bsub -W 08:00 R --file=prostate_Script3.r --save --args 1 550
cmd <- paste(c("bsub -W 01:00 R --file=prostate_Script3.r --save --args ", z, " ", stepsize), collapse="")
system(cmd, intern = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)

waitUntilJobsFinished(1+isMetaRun)


## Script 4
# bsub -W 08:00 R --file=prostate_Script3.r --save --args 1 550
cmd <- paste(c("bsub -W 01:00 R --file=prostate_Script4.r --save"), collapse="")
system(cmd, intern = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)

waitUntilJobsFinished(1+isMetaRun)

delete_successful_lsf(".")
