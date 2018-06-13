#' Drop rows in the SWATH2stats data which have mscores less than or equal to
#' the given mscore threshold.
#'
#' @param data SWATH2stats data structure.
#' @param mscore  Threshold below-which data is unwanted.
#' @param rm.decoy  Drop decoys from the data?
#' @return A copy of the data which has been mscore-filtered.
#' @export
filter_mscore <- function(data, mscore=1.0, rm.decoy=TRUE) {
  if (sum(colnames(data) == "decoy") == 1 & rm.decoy == TRUE) {
    data <- data[data[["decoy"]] == 0, ]
  }

  data.filtered <- data[data[["m_score"]] <= mscore, ]
  message("Original dimension: ", nrow(data), ", new dimension: ", nrow(data.filtered),
          ", difference: ", nrow(data) - nrow(data.filtered), ".")
  return(data.filtered)
}
