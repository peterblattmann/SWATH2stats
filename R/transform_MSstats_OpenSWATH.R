#'  I have no idea what this function does, lets look at it...
#'
#' @param data swath2stats data to modify, I am guessing.
#' @return transformed data?
#' @export
transform_MSstats_OpenSWATH <- function(data) {
  colnames(data) <- gsub("\\.", "", colnames(data))
  colnames(data) <- tolower(colnames(data))

  missing_columns <- setdiff(c("proteinname", "peptidesequence", "precursorcharge",
                               "fragmention", "productcharge", "isotopelabeltype",
                               "condition", "bioreplicate", "run", "intensity"), colnames(data))

  if (length(missing_columns > 0)) {
    stop("The data frame doesn't contain all required columns. The following columns are missing:",
         paste(missing_columns, collapse=", "))
  }

  query_columns <- c("proteinname", "peptidesequence", "precursorcharge", "fragmention",
                     "productcharge", "isotopelabeltype", "condition", "bioreplicate",
                     "run", "intensity")
  add.colnames <- colnames(data)[!(colnames(data) %in% query_column)]

  data <- data[, c(query_columns, add.colnames)]
  new_columns <- gsub(pattern="peptidesequence", replacement="fullpeptidename", x=query_columns)
  colnames(data)[1:10] <- new_columns

  if (length(add.colnames) > 0) {
    message("Additional columns present in the data: ", paste(add.colnames, collapse=", "))
  }

  if (!("mscore" %in% colnames(data))) {
    message("No column 'mscore' present in the data. This column is required for functions estimating FDR or filtering the data.")
  }
  return(data)
}
