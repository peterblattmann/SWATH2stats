#' Filter openSWATH output table according to mscore.
#'
#' This function filters the SWATH data according to the m_score value, as well
#' as to the number of occurence in the data (requant) and within a condition
#' (condition).
#'
#' @param data A data frame containing SWATH data.
#' @param mscore Value that defines the mscore threshold according to which the
#'   data will be filtered.
#' @param rm.decoy  Drop decoys from the data?
#' @return A copy of the data which has been mscore-filtered.
#' @author Peter Blattmann
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore(data, 0.01)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01, 0.8)
#'  data.filtered <- filter_mscore_condition(data, 0.01, 3)
#' }
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
