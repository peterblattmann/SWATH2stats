filter_all_peptides <- function(data){
  
  data.proteins <- gsub("^[[:digit:]]\\/", "", data$ProteinName)
  data.proteins <- unlist(strsplit(as.character(data.proteins), "\\/"))
  data.proteins <- unique(data.proteins)
  
  message("Number of proteins detected: ", length(data.proteins), "\n",
          "Protein identifiers: ", head(data.proteins))
    
  # Delete non-unique features (unique features have "1/" in front of protein name)
  data$ProteinName <- gsub("^1\\/", "", data$ProteinName)
  
  return(data)
}
