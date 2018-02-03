source('common.R')
suppressMessages(library('data.table'))

#' function that calculate the sequence coverage
#'
#' @param path Path to .sam file
#' @param n The length of the reference genome
#' @param m
#' @param sd
#' @return The sequence coverage of the genome associate with the .sam file
insertsStandardDeviation = function(path, n, m, sd) {


  ## vector of changes
  changes = rep(0,n)
  ## deviations vector
  deviations = rep(0,n)

  ## read sam file. Columns needed are:
  ## -col 2:  flag
  ## -col 4:  position
  ## -col 8:  pnext
  ## -col 9:  tlength
  ## -col 10: sequence
  dat = fread(paste('grep "^[^@]" ', path, ' | cut -d"\t" -f2,4,8,9,10', sep = ""), 
    header=FALSE)

  dataLength = dim(dat)[1]

  for(i in 1:dim(dat)[1]) {
    ## Display completness percentage
    displayPercentage(i, dataLength, "Data parsed ")

    x = as.numeric(dat[[1]][i])  ## flag
    y = as.numeric(dat[[2]][i])  ## position
    k = as.numeric(dat[[3]][i])  ## pnext
    v = as.numeric(dat[[4]][i])  ## tlength
    z = nchar(dat[[5]][i])       ## length(sequence)

    ## check flag does not terminate for 11 and length is >0
    ## 11 -> read paired and read paired in proper pair
    if (((bitwAnd(x,3)) == 3) && (v > 0)) {
      changes[y + 1] = abs(k + z - y)
      if (abs(changes[y + 1] - m) > (sd * 2)) {
        deviations[y + 1] = deviations[y + 1] + 1
        deviations[k + z + 1] = deviations[k + z + 1] - 1
      }
    }
  }

  cat("\n")
  rm(list=c("dat"))

  return(deviations)
}