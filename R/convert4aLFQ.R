utils::globalVariables(c("aggregate", "median"))

convert4aLFQ <- function(data, annotation=TRUE){

  if(annotation==TRUE){
    data <- data[, c('ProteinName', 'PeptideSequence', 'FragmentIon', 'NakedSequence', 'PrecursorCharge', 
                     "Intensity", "Condition", "BioReplicate", "Run")]
    data$run_id <- paste(data$Condition, data$BioReplicate, data$Run, sep="_")
  }
  
  colnames(data) <- gsub("ProteinName", "protein_id", colnames(data))
  colnames(data) <- gsub("PrecursorCharge", "precursor_charge", colnames(data))
  colnames(data) <- gsub("FragmentIon", "transition_id", colnames(data))
  colnames(data) <- gsub("Intensity", "transition_intensity", colnames(data))
  colnames(data) <- gsub("PeptideSequence", "peptide_id", colnames(data))
  colnames(data) <- gsub("NakedSequence", "peptide_sequence", colnames(data))
  
  data$transition_id <- gsub("_run[[:digit:]]*$", "", data$transition_id)
  data$transition_id <- paste(data$peptide_sequence, data$transition_id, sep=" ")
  
  data$concentration <- "?"
  
  data <- data[, c("run_id", 'protein_id', 'peptide_id', 'transition_id', 'peptide_sequence', 'precursor_charge', 
                   "transition_intensity", "concentration")]
  
  data.agg <- aggregate(data[,c("transition_id")], by=list(data$peptide_id, data$run_id), length)
  if(median(data.agg$x) == 1){
    warning("The aLFQ package should only be used with transition-level data.")
  }
  return(data)
}
