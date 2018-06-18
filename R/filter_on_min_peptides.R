#' Filter openSWATH output for proteins that are identified by a minimum of n independent peptides.
#'
#' This function removes entries mapping to proteins that are identified by less
#' than n_peptides.
#' Removing single-hit proteins from an analysis can significantly increase the
#' sensitivity under strict protein fdr criteria, as evaluated by
#' e.g. assess_fdr_overall.
#'
#' @param data Data table that is produced by the openSWATH/iPortal workflow.
#' @param n_peptides Number of minimal number of peptide IDs associated with a
#'   protein ID in order to be kept in the dataset.
#' @param rm.decoy  Option to remove the decoys during filtering.
#' @param column which column to use for filtering?
#' @return Returns the filtered data frame with only peptides that map to
#'   proteins with >= n_peptides peptides.
#' @author Moritz Heusel
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01,0.8)
#'  data.max <- filter_on_max_peptides(data.filtered, 5)
#'  data.min.max <- filter_on_min_peptides(data.max, 3)
#' }
#' @export
filter_on_min_peptides <- function(data, n_peptides=6, rm.decoy=TRUE, column="proteinname") {
  data <- unifyProteinGroupLabels(data)
  if(isTRUE(rm.decoy)){
    data <- removeDecoyProteins(data)
  }

  data.prot.pep <- unique(data[, c(column, "fullpeptidename")])
  data.prot.pep.n <- tapply(data.prot.pep[["fullpeptidename"]], data.prot.pep[[column]], length)

  prot.pep.names <- names(data.prot.pep.n[data.prot.pep.n >= n_peptides & !is.na(data.prot.pep.n)])
  data.filtered <- data[data[[column]] %in% prot.pep.names,]

  message("Before filtering:
  Number of proteins: ", length(unique(data[[column]])), "
  Number of peptides: ", length(unique(data[["fullpeptidename"]])), "

Percentage of peptides removed: ", round((length(unique(data[["fullpeptidename"]])) -
                                          length(unique(data.filtered[["fullpeptidename"]]))) /
                                         length(unique(data[["fullpeptidename"]])) * 100, digits=2),
"%", "

After filtering:
  Number of proteins: ", length(unique(data.filtered[[column]])), "
  Number of peptides: ", length(unique(data.filtered[["fullpeptidename"]])))
  return(data.filtered)
}
