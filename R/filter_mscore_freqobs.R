#' Filter openSWATH output table according to mscore.
#'
#' This function filters the SWATH data according to the m_score value, as well
#' as to the number of occurence in the data (requant) and within a condition
#' (condition)
#'
#' @param data  A data frame containing SWATH data.
#' @param mscore  Value that defines the mscore threshold according to which the
#'   data will be filtered.
#' @param percentage  Percentage in which replicas the transition has to reach
#'   the mscore threshold.
#' @param rm.decoy  Option to remove the decoys during filtering.
#' @return Returns a data frame with the filtered data.
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
filter_mscore_freqobs <- function(data, mscore=0.01, percentage=NULL, rm.decoy=TRUE) {
  data[["peptide_charge"]] <- paste(data[["full_peptide_name"]], data[["charge"]])

  if(sum(colnames(data) == "decoy") == 1 & isTRUE(rm.decoy)) {
    data <- data[data[["decoy"]] == 0,]
  }

  data.filtered <- data[data[["m_score"]] <= mscore,]
  data.filtered <- data.table::data.table(data.filtered)

  data.filtered <- data.filtered[, c("peptide_charge", "aggr_peak_area"), with=FALSE]
  data.table::setkey(data.filtered, peptide_charge)
  data.n <- data.filtered[, .N, by="peptide_charge"]

  if (is.null(percentage)) {
    percentage <- 0
  }

  threshold <- nlevels(factor(data[["align_origfilename"]])) * percentage
  message("Peptides need to have been quantified in more conditions than: ",
          threshold, " in order to pass this percentage-based threshold.")

  peptides.filtered <- data.n[data.n[["N"]] >= threshold]
  peptides.filtered <- data.frame("peptide_charge" = peptides.filtered[["peptide_charge"]])

  message("Fraction of peptides selected: ",
          signif(length(unique(peptides.filtered[["peptide_charge"]]))
                 / length(unique(data[["peptide_charge"]])), digits=2))

  data.filtered <- merge(data, peptides.filtered, by="peptide_charge")
  message("Original dimension: ", nrow(data), ", new dimension: ", nrow(data.filtered),
          ", difference: ", nrow(data) - nrow(data.filtered), ".")
  return(data.filtered)
}
