filter_all_peptides <- function(data, column="proteinname", n=6) {
  data.proteins <- gsub("^[[:digit:]]\\/", "", data[[column]])
  data.proteins <- unlist(strsplit(as.character(data.proteins), "\\/"))
  data.proteins <- unique(data.proteins)
  message("Number of proteins detected: ", length(data.proteins))
  message("First ", n, " protein identifiers: ", paste(head(data.proteins, n=n), collapse=", "))
  ## Delete non-unique features (unique features have "1/" in front of protein name)
  data[[column]] <- gsub("^1\\/", "", data[[column]])
  return(data)
}
