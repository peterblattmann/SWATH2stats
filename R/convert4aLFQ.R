convert4aLFQ <- function(data, annotation=TRUE, check_transitions = TRUE) {

  if (annotation==TRUE) {
    data <- data[, c("proteinname", "peptidesequence", "fragmention", "nakedsequence", "precursorcharge",
                     "intensity", "condition", "bioreplicate", "run")]
    data[["run_id"]] <- paste(data[["condition"]], data[["bioreplicate"]], data[["run"]], sep="_")
  }

  colnames(data) <- gsub("proteinname", "protein_id", colnames(data))
  colnames(data) <- gsub("precursorcharge", "precursor_charge", colnames(data))
  colnames(data) <- gsub("fragmention", "transition_id", colnames(data))
  colnames(data) <- gsub("intensity", "transition_intensity", colnames(data))
  colnames(data) <- gsub("peptidesequence", "peptide_id", colnames(data))
  colnames(data) <- gsub("nakedsequence", "peptide_sequence", colnames(data))

  data[["transition_id"]] <- gsub("_run[[:digit:]]*$", "", data[["transition_id"]])
  data[["transition_id"]] <- paste(data[["peptide_sequence"]], data[["transition_id"]], sep=" ")

  data[["concentration"]] <- "?"

  data <- data[, c("run_id", "protein_id", "peptide_id", "transition_id", "peptide_sequence",
                   "precursor_charge", "transition_intensity", "concentration")]
  #check transitions
  data[["protein_id"]] <- factor(data[["protein_id"]])
  data[["peptide_id"]] <- factor(data[["peptide_id"]])
  data[["transition_id"]] <- factor(data[["transition_id"]])

  if (check_transitions) {
    data.agg <- aggregate(data[, c("transition_id")],
                          by=list(data[["peptide_id"]], data[["run_id"]]), length)
    if (median(data.agg[["x"]]) == 1) {
      warning("The aLFQ package should only be used with transition-level data.
              The data only contains one transition per peptide.")
    }
  }

  # convert back to character vector
  data[["protein_id"]] <- as.character(data[["protein_id"]])
  data[["peptide_id"]] <- as.character(data[["peptide_id"]])
  data[["transition_id"]] <- as.character(data[["transition_id"]])
  return(data)
}
