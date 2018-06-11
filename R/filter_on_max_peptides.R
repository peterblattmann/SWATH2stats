filter_on_max_peptides <- function(data, n_peptides, rm.decoy=TRUE, column="proteinname") {
  data <- unifyProteinGroupLabels(data)
  if(isTRUE(rm.decoy)){
    data <- removeDecoyProteins(data)
  }

  data <- data.table::as.data.table(data)
  if (length(grep(column, colnames(data))) > 0) {
    data.table::setnames(data, column, "protein")
  }
  ##data$PEPTIDE <- paste(data$PeptideSequence, data$PrecursorCharge, sep="_")
  data[["peptide"]] <- data[["fullpeptidename"]]
  data.peptides <- data[, c("protein", "peptide", "intensity"), with=FALSE]
  data.table::setkey(data, protein, peptide)

  data.peptides.int <- data.peptides[, sum(intensity), by="protein,peptide"]
  data.table::setnames(data.peptides.int, "V1", "sum.intensity")

  data.table::setkey(data.peptides.int, protein)
  data.peptides.int <- data.peptides.int[order(data.peptides.int[["sum.intensity"]], decreasing=TRUE), ]
  peptides.sel <- unique(data.peptides.int[, head(.SD, n_peptides), by=protein])
  data.filtered <- data.frame(data[peptide %in% peptides.sel[["peptide"]], ])

  message("Before filtering: ", "\n",
          "  Number of proteins: ", length(unique(data[["protein"]])), "\n",
          "  Number of peptides: ", length(unique(data[["peptide"]])), "\n\n",
          "Percentage of peptides removed: ", round(
          (length(unique(data[["peptide"]])) - length(unique(data.filtered[["peptide"]]))) /
          length(unique(data[["peptide"]])) * 100, digits=2), "%", "\n\n",
          "After filtering: ", "\n",
          "  Number of proteins: ", length(unique(data.filtered[["protein"]])), "\n",
          "  Number of peptides: ", length(unique(data.filtered[["peptide"]])))

  colnames(data.filtered) <- gsub("protein", "proteinname", colnames(data.filtered))
  data.filtered <- data.filtered[, -which(colnames(data.filtered) == "peptide")]
  return(data.filtered)
}
