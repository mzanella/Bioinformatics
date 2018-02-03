suppressMessages(library('data.table'))

#' function that calculate the length of genome inserts
#'
#' @param path Path to .sam file
#' @param n The length of the reference genome
#' @return The length of genome inserts
genomeInserts = function(path, n) {

  ## vector of changes
  changes = rep(0,n)

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

    x = as.numeric(dat[[1]][i]) ## flag
    y = as.numeric(dat[[2]][i]) ## position
    v = as.numeric(dat[[4]][i]) ## tlegth
    ## pnext + legth(sequence)
    z = as.numeric(dat[[3]][i]) + nchar(dat[[5]][i])

    ## check flag does not terminate for 11 and length is >0
    ## 11 -> read paired and read paired in proper pair
    if (((bitwAnd(x,3)) == 3) && (v > 0)) {
        changes[y + 1] = abs(z - y)
    }
  }

  cat("\n")
  rm(list=c("dat"))

  return(changes)
}