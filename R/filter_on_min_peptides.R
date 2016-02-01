filter_on_min_peptides <- function(data, n_peptides){
  data.prot.pep <- unique(data[,c("ProteinName", "FullPeptideName")])

  data.prot.pep.n <- tapply(data.prot.pep$FullPeptideName, data.prot.pep$ProteinName, length)

  prot.pep.names <- names(data.prot.pep.n[data.prot.pep.n >= n_peptides & !is.na(data.prot.pep.n)])
  data.filtered <- data[data$ProteinName %in% prot.pep.names,]

  message("Before filtering: ", "\n",
          "  Number of proteins: ", length(unique(data$ProteinName)), "\n",
          "  Number of peptides: ", length(unique(data$FullPeptideName)), "\n\n",
          "Percentage of peptides removed: ", round((length(unique(data$FullPeptideName)) - length(unique(data.filtered$FullPeptideName)))/length(unique(data$FullPeptideName))*100, digits=2), "%", "\n\n",
          "After filtering: ", "\n",
          "  Number of proteins: ", length(unique(data.filtered$ProteinName)), "\n",
          "  Number of peptides: ", length(unique(data.filtered$FullPeptideName)), "\n")


  return(data.filtered)
}