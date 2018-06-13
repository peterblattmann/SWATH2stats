#' Bring together protein group labels.
#'
#' @param data  SWATH2stats data to modify.
#' @param column Which column to use for unifying the groups.
#' @return slightly modified data structure.
#' @export
unifyProteinGroupLabels <- function(data, column="proteinname") {
  data[[column]] <- as.character(data[[column]])
  ids <- grep("^([2-9])|([1-9][0-9][0-9]*)/", data[[column]])
  identifiers <- data[ids, column]
  identifiers_split <- strsplit(as.character(identifiers), "/")
  identifiers_split_sorted <- lapply(identifiers_split, function(x) { sort(x) })
  identifiers_sorted <- sapply(identifiers_split_sorted, function(x) { paste(x, collapse="/") } )
  data[ids, column] <- identifiers_sorted
  return(data)
}
