unifyProteinGroupLabels <- function(data){
  data$ProteinName <- as.character(data$ProteinName)
  ids <- grep("^([2-9])|([1-9][0-9][0-9]*)/", data$ProteinName)
  identifiers <- data[ids,"ProteinName"]
  identifiers_split <- strsplit(as.character(identifiers), "/")
  identifiers_split_sorted <- lapply(identifiers_split, function(x){sort(x)})
  identifiers_sorted <- sapply(identifiers_split_sorted, function(x){paste(x, collapse="/")})
  data[ids, "ProteinName"] <- identifiers_sorted
  return(data)
}
