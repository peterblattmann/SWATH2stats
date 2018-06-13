#' Write out a matrix of peptides from swath2stats data.
#'
#' @param data SWATH2stats data to write out.
#' @param write.csv  Write the result to a csv file?  This should be folded into the next arg.
#' @param fun.aggregate  What function to use when aggregating the set of intensities?
#' @param filename  Where to write the data?
#' @param rm.decoy  Remove the decoys?
#' @return the peptides as a matrix!
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
