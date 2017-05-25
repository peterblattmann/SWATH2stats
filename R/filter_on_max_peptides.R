utils::globalVariables(c("PROTEIN", "PEPTIDE", "Intensity", ".SD", "head"))

filter_on_max_peptides <- function(data, n_peptides, rm.decoy = TRUE){
  data <- unifyProteinGroupLabels(data)
  if(isTRUE(rm.decoy)){
    data <- removeDecoyProteins(data)
  }
  
  data<- data.table(data)
  if(length(grep("ProteinName", colnames(data))) > 0){
    setnames(data, "ProteinName", "PROTEIN")
  }
  
  #data$PEPTIDE <- paste(data$PeptideSequence, data$PrecursorCharge, sep="_")
  data$PEPTIDE <- data$FullPeptideName


  data.peptides <- data[,c("PROTEIN", "PEPTIDE", "Intensity"), with=FALSE]
  setkey(data, PROTEIN, PEPTIDE) 

  data.peptides.int <- data.peptides[, sum(Intensity), by="PROTEIN,PEPTIDE"]
  setnames(data.peptides.int, "V1", "SUM.INTENSITY")

  setkey(data.peptides.int, PROTEIN)
  data.peptides.int <- data.peptides.int[order(data.peptides.int$SUM.INTENSITY, decreasing=TRUE), ]

  peptides.sel <- unique(data.peptides.int[, head(.SD, n_peptides), by=PROTEIN])

  data.filtered <- data.frame(data[PEPTIDE %in% peptides.sel$PEPTIDE,])


  message("Before filtering: ", "\n",
          "  Number of proteins: ", length(unique(data$PROTEIN)), "\n",
          "  Number of peptides: ", length(unique(data$PEPTIDE)), "\n\n",
          "Percentage of peptides removed: ", round((length(unique(data$PEPTIDE)) - length(unique(data.filtered$PEPTIDE)))/length(unique(data$PEPTIDE))*100, digits=2), "%", "\n\n",
          "After filtering: ", "\n", 
          "  Number of proteins: ", length(unique(data.filtered$PROTEIN)), "\n",
          "  Number of peptides: ", length(unique(data.filtered$PEPTIDE)), "\n")

  colnames(data.filtered) <- gsub("PROTEIN", "ProteinName", colnames(data.filtered))
  data.filtered <- data.filtered[,-which(colnames(data.filtered) == "PEPTIDE")]

  return(data.filtered)
}


