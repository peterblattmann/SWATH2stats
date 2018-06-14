#' Filter data based on a given FDR and mscore. (I think)
#'
#' @param data  SWATH2stats data structure to filter.
#' @param mscore  Chosen mscore limit.
#' @param percentage  Look for peptides in > this percent of the possible
#'   conditions, if they are found in less than this threshold, bye bye!
#' @param rm.decoy  Drop the decoys!?
#' @return filtered data!
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
