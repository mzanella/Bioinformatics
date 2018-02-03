suppressMessages(library('data.table'))

#' Shows the help for the script
#'
showHelp = function() {
  message("
Script created for the Bioinformatics project

Usage: 

Rscript main.R genome_length path command_type [mean] [standard_deviation]
- genome_length: number, the lenght of the genome
- path: relative path to the sam file
- command_type: it could be
    * 'physical' - calculate the physical coverage
    * 'sequence' - calculate the sequence coverage
    * 'insert' - calculate the length of the inserts
    * 'deviation' - calculate the standard deviation
- mean: needed only in order to calculate the 'deviation'
- standard_deviation: needed only in order to calculate the 'deviation'
          
Rscript main.R [-h|--help]
Shows this help
")
}


#' Checks if the minimun number of arguments is given
#'
checkMinArgs = function() {
  if (length(args)<3) {
    stop(paste("At least 3 argument must be supplied. Run ",
      "'Rscript main.R -h' for more details", sep=""), call.=FALSE)
  }
}


#' Checks if the command typed from the user is correct
#'
checkCommand = function() {
  if (!(command %in% c("physical", "sequence", "deviation", "insert", "all"))) {
    stop(paste("The command must be 'physical', 'sequence', 'insert', ",
      "'deviation' or 'all'. Run 'Rscript main.R -h' for more details", sep=""),
    call.=FALSE)
  }
}


#' Generate data to be saved in the .wig file
#'
#' @param command The type of data to save
#' @param data The data calculated
#' @param n The length of the genome
#' @return Data formatted to be saved in the .wig file
generateWigData = function(command, data, n) {
  head = "fixedStep chrom=genome start=1 step=1 span=1\n\n"
  toPrint = ""

  if (command %in% c('physical', 'sequence', 'deviation')) {
    toPrint = paste(cumsum(data), collapse="\n")
    toPrint = paste(head, toPrint, sep="")
  } else if (command == 'insert') {
    toPrint = paste(data, collapse="\n")
    toPrint = paste(head, toPrint, sep="")
  }

  return(toPrint)
}


#' Save data to file
#'
#' @param filename Name of the file
#' @param data Data to save
saveFile = function(filename, data) {
  library(data.table)
  fwrite(list(data), filename, sep=" ", quote=FALSE)
}


#' Connection to a file that must be read
#'
#' @param path Path to the file
#' @return The connection to the file
readFile = function(path) {
  return(file(description=path, open="r"))
}


#' Splitting by tab character
#'
#' @param string String that must be splitted
#' @return List of string retrieved splitting the argument given
splitByTabs = function(string) {
  return(strsplit(string, "\t")[[1]])
}


#' Calculate the mean discarding 0s
#'
#' @param dat Data of which calculate the mean value
#' @return The mean value
meanValue = function(dat) {
  return(mean(removeValue(dat,0)))
}


#' Calculate the standard deviation discarding 0s
#'
#' @param dat Data of which calculate the standard deviation
#' @return The standard deviation
standardDeviation = function(dat) {
  return(sd(removeValue(dat, 0)))
}


#' Funtion that removes all field equals to a certain value
#'
#' @param dat Data to which remove a certain value
#' @param toRemove Value to remove
removeValue = function(dat, toRemove = 0) {
  return(dat[dat != toRemove])
}


#' Function to read single column data file efficiently
#'
#' @param path Path to the file to read
#' @param comment Char that represent the comment, default is the hashtag symbol
#' @return Array of data read
readSingleColumnFile = function(path, comment = "#") {
  dat = fread(paste('grep "^[^', comment, ']" ', path, ' | cut -d"\t" -f1', 
    sep = ""), header=FALSE)
  dat = sapply(dat, as.numeric)
  return(dat)
}

#' Function to display the percentage of completeness updating the same row
#'
#' @param actual Actual value completed
#' @param total Value to reach
#' @param string String to display
displayPercentage = function(actual, total, string = "") {
  command = paste("\r", string, sep = "")
  cat(command, format(round((actual*100)/total, 2), nsmall = 2), "%")
}