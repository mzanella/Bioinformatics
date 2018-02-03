source("common.R")
source("physicalcoverage.R")
source("sequencecoverage.R")
source("genomeinserts.R")
source("genomestandarddeviation.R")

allValues = function(path, n) {

  physical = rep(0,n)
  sequence = rep(0,n)
  inserts = rep(0,n)

  ## connection to the file in read-only mode
  f = file(description=path, open="r")
  
  while (length(line <- readLines(f, n = 1, warn = FALSE)) > 0) {
    ## if it is not a comment line
    if (!(startsWith(line, '@'))) {
      cols = strsplit(line, "\t")[[1]]
      ## check if alignments are correct and their length is >0
      physical = physicalCalc(cols, physical)
      sequence = sequenceCalc(cols, sequence)
      inserts = insertsCalc(cols, inserts)
    }
  }

  close(f)

  saveFile("physical.wig", generateWigData("physical", physical, n))
  saveFile("sequence.wig", generateWigData("sequence", sequence, n))
  saveFile("insert.wig", generateWigData("insert", inserts, n))

  rm(list=c("physical", "sequence", "inserts"))

  dat = readSingleColumnFile('insert.wig', 'f')
  saveFile("deviation.wig", 
    generateWigData('deviation', 
      insertsStandardDeviation(path, n, meanValue(dat), standardDeviation(dat))
      , n)
  )

  return(changes)
}

physicalCalc = function(cols, changes) {
  if (((bitwAnd(as.numeric(cols[2]),3)) == 3) && (as.numeric(cols[9]) > 0)) {
    changes[as.numeric(cols[4]) + 1] = changes[as.numeric(cols[4]) + 1] + 1
    changes[as.numeric(cols[8]) + nchar(cols[10]) + 1] = 
      changes[as.numeric(cols[8]) + nchar(cols[10]) + 1] - 1
  }
  return(changes)
}

sequenceCalc = function(cols, changes) {
  if ((bitwAnd(as.numeric(cols[2]),5)) != 5) {
    changes[as.numeric(cols[4]) + 1] = changes[as.numeric(cols[4]) + 1] + 1
    changes[as.numeric(cols[4]) + nchar(cols[10]) + 1] = 
      changes[as.numeric(cols[4]) + nchar(cols[10]) + 1] - 1
  }
  return(changes)
}

insertsCalc = function(cols, changes) {
  if (((bitwAnd(as.numeric(cols[2]),3)) == 3) && (as.numeric(cols[9]) > 0)) {
    changes[as.numeric(cols[4]) + 1] = 
      abs(as.numeric(cols[8]) + nchar(cols[10]) - as.numeric(cols[4]))
  }
  return(changes)
}