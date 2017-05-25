utils::globalVariables(c("head"))

filter_proteotypic_peptides <- function(data, rm.decoy = TRUE){
  
  data <- unifyProteinGroupLabels(data)
  if(isTRUE(rm.decoy)){
    data <- removeDecoyProteins(data)
  }
  
  data.proteins <- gsub("^[[:digit:]]\\/", "", data$ProteinName)
  data.proteins <- unlist(strsplit(as.character(data.proteins), "\\/"))
  data.proteins <- unique(data.proteins)
  
  message("Number of proteins detected: ", length(data.proteins), "\n",
          "Protein identifiers: ", paste(head(data.proteins), collapse=", "))
    
  # Delete non-unique features (unique features have "1/" in front of protein name)
  data.unique <- data[grep("^1/", data$ProteinName), ]
  data.unique$ProteinName <- gsub("^1\\/", "", data.unique$ProteinName)
  message("Number of proteins detected that are supported by a proteotypic peptide: ", 
          length(unique(data.unique$ProteinName)), "\n",
          "Number of proteotypic peptides detected: ", 
          length(unique(data.unique$FullPeptideName)))
  
  return(data.unique)
}