#' Filter only for the highest intense peptides
#'
#' In order to reduce the data, the data is filtered only for the proteins with
#' the highest intensity peptides.
#'
#' @param data A data frame containing SWATH data with the column names:
#'   ProteinNames, PeptideSequence, PrecursorCharge, Intensity.
#' @param n_peptides Maximum number of highest intense peptides to filter the
#'   data on.
#' @param rm.decoy  Option to remove the decoys during filtering.
#' @param column which column to use for filtering?
#' @return  Returns a data frame of the filtered data.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01,0.8)
#'  data.max <- filter_on_max_peptides(data.filtered, 5)
#' @export
filter_on_max_peptides <- function(data, n_peptides=6, rm.decoy=TRUE, column="proteinname") {
  data <- unifyProteinGroupLabels(data)
  if (isTRUE(rm.decoy)) {
    data <- removeDecoyProteins(data)
  }

  data <- data.table::as.data.table(data)
  if (length(grep(column, colnames(data))) > 0) {
    data.table::setnames(data, column, "protein")
  }
  ##data$PEPTIDE <- paste(data$PeptideSequence, data$PrecursorCharge, sep="_")
  ## I have absolutely no clue why, but if this is left as its default factor, then
  ## when we do the setkey(data, protein, peptide) in a couple lines, the
  ## peptide names will be messed up.  To me, this feels like a bug in R, as
  ## there is no reason for it to reorder the elements in the factor when doing
  ## the setkey(), yet it clearly is.
  data[["peptide"]] <- as.character(data[["fullpeptidename"]])

  data.peptides <- data[, c("protein", "peptide", "intensity"), with=FALSE]
  ## This line is where my copy of data changes w.r.t. the original code.
  ##data.table::setkey(data, protein, peptide)
  data.table::setkeyv(x=data, cols=c("protein", "peptide"))

  ## I am not sure how to get around this bit of odd NSE.
  intensity <- NULL
  data.peptides.int <- data.peptides[, sum(intensity), by="protein,peptide"]
  data.table::setnames(data.peptides.int, "V1", "sum.intensity")

  ##data.table::setkey(data.peptides.int, protein)
  data.table::setkeyv(x=data.peptides.int, cols="protein")
  data.peptides.int <- data.peptides.int[order(data.peptides.int[["sum.intensity"]], decreasing=TRUE), ]

  ## .SD is data.table shorthand for _s_ubset _d_atatable.
  .SD <- NULL
  peptides.sel <- unique(data.peptides.int[, head(x=.SD, n=n_peptides), by="protein"])
  selected_idx <- data[["peptide"]] %in% peptides.sel[["peptide"]]
  ##data.filtered <- data.frame(data[peptide %in% peptides.sel[["peptide"]], ])
  data.filtered <- as.data.frame(data[selected_idx, ])

  message("Before filtering:\n",
          "  Number of proteins: ", length(unique(data[["protein"]])), "\n",
          "  Number of peptides: ", length(unique(data[["peptide"]])), "\n\n",
          "Percentage of peptides removed: ", round(
          (length(unique(data[["peptide"]])) - length(unique(data.filtered[["peptide"]]))) /
          length(unique(data[["peptide"]])) * 100, digits=2), "%\n\n",
          "After filtering:\n",
          "  Number of proteins: ", length(unique(data.filtered[["protein"]])), "\n",
          "  Number of peptides: ", length(unique(data.filtered[["peptide"]])))

  colnames(data.filtered) <- gsub("protein", "proteinname", colnames(data.filtered))
  data.filtered <- data.filtered[, -which(colnames(data.filtered) == "peptide")]
  return(data.filtered)
}
