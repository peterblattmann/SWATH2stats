removeDecoyProteins <- function(data, column="proteinname") {
  ids <- grep("^[1-9][0-9]*[0-9]*/.*DECOY.*", data[[column]])
  identifiers <- data[ids, column]
  identifiers_split <- strsplit(as.character(identifiers), "/")
  identifiers_split_removed <- lapply(identifiers_split, rmDecoyProt)
  identifiers_removed <- sapply(identifiers_split_removed, function(x) { paste(x, collapse="/") })
  data[ids, column] <- identifiers_removed
  return(data)
}

rmDecoyProt <- function(x, pattern="DECOY") {
  ids <- grep(pattern=pattern, x=x)
  ids.sel <- grep(pattern=pattern, x=x, invert=TRUE)
  x[1] <- as.numeric(x[1]) - length(ids)
  x <- x[ids.sel]
  return(x)
}
