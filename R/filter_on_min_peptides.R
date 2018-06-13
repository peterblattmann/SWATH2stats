#' Filter-out too-low peptides from the data!
#'
#' @param data SWATH2stats data to filter.
#' @param n_peptides Number of first n peptides to keep
#' @param rm.decoy  Drop the decoys?
#' @param column which column to use for filtering?
#' @return filtered data!
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
