source("common.R")
source("physicalcoverage.R")
source("sequencecoverage.R")
source("genomeinserts.R")
source("genomestandarddeviation.R")
#source("all.R")

suppressMessages(require("compiler"))
suppressMessages(compiler::enableJIT(3))

args=commandArgs(trailingOnly=TRUE)
if (args[1] %in% c("-h", "--help")) {
  showHelp()
} else {

  checkMinArgs()

  n = args[1]
  path = args[2]
  command = tolower(args[3])
  m = -1
  sd = -1

  checkCommand()

  if (length(args)<5 && command=="deviation") {
    warning(paste("'deviation' command need 5 argument, but less are given.",
      " Mean and standard deviation will be calculate from 'insert.wig' file",
      " if it is present", sep=""), call.=FALSE)
  } else if (toupper(command)=="deviation") {
    m = args[4]
    sd = args[5]
  }

  res=""

  switch(command,
    physical = {
      res = physicalCoverage(path, n)
     },
    sequence = {
      res = sequenceCoverage(path, n)
    },
    deviation = {
      if (m == -1 && sd == -1) {
        dat = readSingleColumnFile('insert.wig', 'a-zA-Z')
        m = round(meanValue(dat), digits = 0)
        sd = standardDeviation(dat)
        print(paste("Mean:", m))
        print(paste("Standard deviation:", sd))
      }
      res = insertsStandardDeviation(path, n, m, sd)
    },
    insert = {
      res = genomeInserts(path, n)
    }
  )

  if(!(command %in% c("all"))) {
    toPrint = generateWigData(command, res, n)
    saveFile(paste(command, ".wig", sep=""), toPrint)
  }
}