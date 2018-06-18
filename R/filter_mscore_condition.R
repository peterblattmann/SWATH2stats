#' Filter openSWATH output table according to mscore.
#'
#' This function filters the SWATH data according to the m_score value, as well
#' as to the number of occurence in the data (requant) and within a condition
#' (condition).
#'
#' @param data A data frame containing SWATH data.
#' @param mscore Value that defines the mscore threshold according to which the
#'   data will be filtered.
#' @param n.replica Number of measurements within at least one condition that
#'   have to pass the mscore threshold for this transition.
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
filter_mscore_condition <- function(data, mscore=1.0, n.replica, rm.decoy=TRUE) {

  ## a column with unique identifiers for each Precursor (e.g. ADFSDF 3) is generated
  data[["peptide_charge"]] <- paste(data[["fullpeptidename"]], data[["charge"]])
  ## a column with unique identifiers for each Precursor and Condition is genrated ADFSDF 3 Condition 1
  data[["peptide_charge_condition"]] <- paste(data[["peptide_charge"]], data[["condition"]])

  # decoys are removed if present
  if(sum(colnames(data) == "decoy") == 1 & isTRUE(rm.decoy)){
    data <- data[data[["decoy"]] == 0,]
    #data <- subset(data, decoy == 0)
  }

  # only data that is below the indicated mscore is selected and then only the unique data selected
  #data.filtered <- subset(data, m_score <= mscore)
  data.filtered <- data[data[["m_score"]] <= mscore,]
  data.filtered <- unique(data.filtered[,c("peptide_charge", "peptide_charge_condition", "aggr_peak_area")])
  data.filtered <- data.table::data.table(data.filtered)

  data.table::setkey(data.filtered, peptide_charge, peptide_charge_condition, aggr_peak_area)
  # number of occurences of Precursor per Condition is calculated
  data.n <- data.filtered[, .N, by="peptide_charge,peptide_charge_condition"]

  # only precursors that are present in more that n.replica are selected
  precursor.filtered <- data.n[data.n[["N"]] >= n.replica]
  precursor.filtered <- data.frame(peptide_charge = unique(precursor.filtered[["peptide_charge"]]))

  message("Fraction of peptides selected: ",
          signif(length(unique(precursor.filtered[["peptide_charge"]])) /
                 length(unique(data[["peptide_charge"]])), digits=2))

  # only data that is present in the precursor.filtered list is selected
  data.filtered <- merge(data, precursor.filtered, by.x="peptide_charge", by.y="peptide_charge")
  message("Original dimension: ", nrow(data), ", new dimension: ", nrow(data.filtered),
          ", difference: ", nrow(data) - nrow(data.filtered), ".")
  return(data.filtered)
}
