write_matrix_peptides <- function(data, write.csv=FALSE,
                                  filename="SWATH2stats_overview_matrix_peptidelevel.csv",
                                  rm.decoy=FALSE) {
  if (rm.decoy == TRUE) {
    data <- subset(data, data[["decoy"]] == 0)
  }
  data.peptide <- data[, c("proteinname", "run_id", "fullpeptidename", "intensity")]
  ProteinName_FullPeptideName <- paste(data.peptide[["proteinname"]],
                                       data.peptide[["fullpeptidename"]], sep="_")
  data.peptide <- cbind(ProteinName_FullPeptideName, data.peptide)
  data.peptide.table <- dcast(data.peptide, ProteinName_FullPeptideName ~ run_id,
                              value.var="intensity", fun.aggregate=sum)
  if (isTRUE(write.csv)) {
    write.csv(data.peptide.table, file=filename, row.names=FALSE, quote=FALSE)
    message("Peptide overview matrix ", filename," written to working folder.")
  }
  return(data.peptide.table)
}
