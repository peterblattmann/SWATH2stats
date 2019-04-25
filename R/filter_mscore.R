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
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore(data, 0.01)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01, 0.8)
#'  data.filtered <- filter_mscore_condition(data, 0.01, 3)
#' @export
filter_mscore <- function(data, mscore, rm.decoy = TRUE, mscore.col = "m_score") {
    if (sum(colnames(data) == "decoy") == 1 & rm.decoy == TRUE) {
        data <- data[data$decoy == 0, ]
        # subset(data, decoy == 0)
    }
    mscore.col <- JPP_update(data, mscore.col)
    data.filtered <- data[data[, mscore.col] <= mscore, ]
    message("Dimension difference: ", paste(dim(data) - dim(data.filtered), collapse = ", "))
    return(data.filtered)
}
