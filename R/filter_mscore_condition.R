#' Filter SWATH2stats data based on mscore/condition
#'
#' @param data SWATH2stats data.
#' @param mscore Threshold below which to drop rows.
#' @param n.replica  Minimum required number of replicates.
#' @param rm.decoy  Drop decoy rows?
#' @return filtered data structure
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
