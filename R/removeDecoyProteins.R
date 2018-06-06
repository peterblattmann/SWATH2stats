removeDecoyProteins <- function(data){
  ids <- grep("^[1-9][0-9]*[0-9]*/.*DECOY.*", data$ProteinName)
  identifiers <- data[ids,"ProteinName"]
  identifiers_split <- strsplit(as.character(identifiers), "/")
  identifiers_split_removed <- lapply(identifiers_split, rmDecoyProt)
  identifiers_removed <- vapply(identifiers_split_removed, function(x){paste(x, collapse="/")}, "a")
  data[ids,"ProteinName"] <- identifiers_removed
  return(data)
}

rmDecoyProt <- function(x){
  ids <- grep("DECOY", x)
  ids.sel <- grep("DECOY", x, invert=TRUE)
  x[1] <- as.numeric(x[1])-length(ids)
  x <- x[ids.sel]
  return(x)
}
