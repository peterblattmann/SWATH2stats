#' This function is named filter_all_peptides, but it appears to me that it
#' primarily serves to remove the leading slash number characters from the
#' peptide names -- which is certainly a good thing(TM), but is that filtering?
#' Perhaps I am still missing something
#'
#' @param data  SWATH2stats data structure to be modified.
#' @param column  Which column contains the data to modify?
#' @param n  How many new IDs should we print to show if this worked?
#' @return  Slightly modified data in the same format.
#' @export
filter_all_peptides <- function(data, column="proteinname", n=6) {
  data_proteins <- gsub(pattern="^[[:digit:]]\\/", replacement="", x=data[[column]])
  data_proteins <- unlist(strsplit(as.character(data_proteins), "\\/"))
  data_proteins <- unique(data_proteins)
  message("Number of proteins detected: ", length(data_proteins))
  message("First ", n, " protein identifiers: ", paste(head(data_proteins, n=n), collapse=", "))
  ## Delete non-unique features (unique features have "1/" in front of protein name)
  data[[column]] <- gsub(pattern="^1\\/", replacement="", x=data[[column]])
  return(data)
}
