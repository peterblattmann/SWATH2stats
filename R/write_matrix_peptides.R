#' Writes out an overview matrix of peptides mapping to a FDR quality controlled
#' protein master list at controlled global peptide FDR quality.
#'
#' Writes out an overview matrix on peptide level of a supplied (unfiltered or
#' prefiltered) OpenSWATH results data frame.
#' The peptide quantification is achieved by summing the areas under all 6
#' transitions per precursor and summing all precursors per FullPeptideName.
#' In order to keep the peptide-to-protein association, the FullPeptideName is
#' joined with the ProteinName.
#'
#' @param data A data frame containing annotated OpenSWATH/pyProphet data.
#' @param write.csv Option to determine if table should be written automatically
#'   into csv file.
#' @param fun.aggregate  What function to use when aggregating the set of
#'   intensities?
#' @param filename File base name of the .csv matrix written out to the working
#'   folder.
#' @param rm.decoy Logical whether decoys will be removed from the data
#'   matrix. Defaults to FALSE. It's sometimes useful to know how decoys behave
#'   across a dataset and how many you allow into your final table with the
#'   current filtering strategy.
#' @return the peptides as a matrix! also output .csv matrix is written to the
#'   working folder.
#' @author Moritz Heusel
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  write_matrix_peptides(data)
#' }
#' @export
write_matrix_peptides <- function(data, write.csv=FALSE,
                                  fun.aggregate=sum,
                                  filename="SWATH2stats_overview_matrix_peptidelevel.csv",
                                  rm.decoy=FALSE) {
  if (rm.decoy == TRUE) {
    data <- subset(data, data[["decoy"]] == 0)
  }
  data.peptide <- data[, c("proteinname", "run_id", "fullpeptidename", "intensity")]
  ProteinName_FullPeptideName <- paste(data.peptide[["proteinname"]],
                                       data.peptide[["fullpeptidename"]], sep="_")
  data.peptide <- cbind(ProteinName_FullPeptideName, data.peptide)
  data.peptide.table <- reshape2::dcast(data.peptide, ProteinName_FullPeptideName ~ run_id,
                                        value.var="intensity", fun.aggregate=fun.aggregate)
  if (isTRUE(write.csv)) {
    write.csv(data.peptide.table, file=filename, row.names=FALSE, quote=FALSE)
    message("Peptide overview matrix ", filename, " written to working folder.")
  }
  return(data.peptide.table)
}
