#' Convert table into the format for aLFQ
#'
#' This functions selects the columns necessary for the aLFQ R package.
#'
#' @param data A data frame containing the SWATH data in transition-level format
#' @param annotation Option to indicate if the data has been annotated, i.e. if
#'   the columns Condition, Replicate, Run are present. If option is set to true
#'   it will write a new run_id as a string of the combination of these three
#'   columns.
#' @param  check_transitions if number of transitions should be checked. As input only
#'   transition-level data should be used and therefore this is
#'   checked. However, this makes the function slow and herewith be omitted.
#' @return Returns a data frame in the appropriate format for aLFQ.
#' @references Rosenberger G, Ludwig C, Rost HL, Aebersold R, Malmstrom L. aLFQ:
#'   an R-package for estimating absolute protein quantities from label-free
#'   LC-MS/MS proteomics data. Bioinformatics. 2014 Sep 1;30(17):2511-3. doi:
#'   10.1093/bioinformatics/btu200.
#' @author Peter Blattmann
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#'  data.aLFQ <- convert_aLFQ(raw)
#' }
#' @export
convert_aLFQ <- function(data, annotation=TRUE, check_transitions=TRUE) {

  if (isTRUE(annotation)) {
    data <- data[, c("proteinname", "peptidesequence", "fragmention",
                     "nakedsequence", "precursorcharge", "intensity",
                     "condition", "bioreplicate", "run")]
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
