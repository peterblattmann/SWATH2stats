#' Count the number of analytes observed in a swath2stats data structure.
#'
#' @param data  swath2stats data structure.
#' @param column_levels  Use these columns to collate the data with stats::aggregate.
#' @param column_by Use this column for collecting the resulting counts.
#' @param rm_decoy  Keep the decoys?
#' @return  data frame of counted peptides.
#' @export
count_analytes <- function(data,
                           column_levels=c("transition_group_id", "fullpeptidename", "proteinname"),
                           column_by="run_id", rm_decoy=TRUE) {
  if (sum(colnames(data) == "decoy") == 1 & isTRUE(rm_decoy)) {
    data <- data[data[["decoy"]] == 0, ]
  }

  data_n <- aggregate(data[, column_levels],
                      by=list(data[, column_by]),
                      function(x) length(unique(x)))
  colnames(data_n)[1] <- column_by
  return(data_n)
}
