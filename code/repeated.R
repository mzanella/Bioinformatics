require(data.table)

#' Save on `more_than_2_alignments.csv` and `only_1_alignments` files reads
#' found, respectively 3 or more times and only twice
#' 
#' @param path Path to sam file create with bwa mem -a
findRepeatedReads = function(path) {
  ## takes only the identifier 

  dat = fread(paste('grep "^[^@]" ', path ,' | cut -d"\t" -f1', sep = ""),
    header=FALSE)[[1]]

  ## create a table with the identifier -> calculate the frequency of each row
  t = table(dat)

  ## extract identifier found 3 or more times
  m2 = t[t>2]

  ## extract identifier found twice
  m1 = t[t==2]

  ## write results on files
  write.table(m2, file='more_than_2_alignements.csv', col.names=FALSE,
    row.names=FALSE)
  write.table(m2, file='only_1_alignement.csv', col.names=FALSE, row.names=FALSE)
}

## catch the user input and call the function
args=commandArgs(trailingOnly=TRUE)
findRepeatedReads(args[1])