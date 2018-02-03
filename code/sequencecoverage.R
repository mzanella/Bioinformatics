source('common.R')
suppressMessages(library('data.table'))
#' function that calculate the sequence coverage
#'
#' @param path Path to .sam file
#' @param n The length of the reference genome
#' @return The sequence coverage of the genome associate with the .sam file
sequenceCoverage = function(path, n) {

  ## vector of changes
  changes = rep(0,n)
  
  ## read sam file. Columns needed are:
  ## -col 2:  flag
  ## -col 4:  position
  ## -col 10: sequence
  dat = fread(paste('grep "^[^@]" ', path, ' | cut -d"\t" -f2,4,10', sep = ""), 
    header=FALSE)

  dataLength = dim(dat)[1]

  for(i in 1:dim(dat)[1]) {
    ## Display completness percentage
    displayPercentage(i, dataLength, "Data parsed ")

    ## flag
    x = as.numeric(dat[[1]][i])
    ## position + 1
    y = as.numeric(dat[[2]][i]) + 1
    ## position + 1 + length(sequence)
    z = y + nchar(dat[[3]][i])

    ## check flag does not terminate for 101
    ## 101 -> read paired and read unmapped
    if ((bitwAnd(x,5)) != 5) {
      changes[y] = changes[y] + 1
      changes[z] = changes[z] - 1
    }
  }

  cat("\n")
  rm(list=c("dat"))

  return(changes)
}