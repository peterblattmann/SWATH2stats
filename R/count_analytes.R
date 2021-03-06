#' Counts the analytes across the different injections
#'
#' This functions counts the number of different peakgroups, peptides and
#' proteins in different injections.
#'
#' @param data A data frame containing SWATH data.
#' @param column_levels  Columns in which the number of unique identifiers 
#'   should be counted.
#' @param column_by Column for which the different identifiers should be counted
#'   for, e.g. for the different injections.
#' @param rm_decoy Option to not remove decoy before counting.
#' @return Returns a data frame with the count of the different identifiers per
#'   e.g. injection.
#' @author Peter Blattmann
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  counts <- count_analytes(data)
#'  }
#' @importFrom stats aggregate
#' @export
count_analytes <- function(data, 
                           column_levels = c("transition_group_id", 
                                             "FullPeptideName", 
                                             "ProteinName"), 
                           column_by = "run_id", 
                           rm_decoy = TRUE) {
    if (sum(colnames(data) == "decoy") == 1 & isTRUE(rm_decoy)) {
        # data <- subset(data, decoy == 0)
        data <- data[data$decoy == 0, ]
    }

    data.n <- aggregate(data[, column_levels], 
                        by = list(data[, column_by]), 
                        function(x) length(unique(x)))
    colnames(data.n)[1] <- column_by
    return(data.n)
}
