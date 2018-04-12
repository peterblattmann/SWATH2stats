convert4pythonscript <- function(data, replace.Unimod = TRUE){

  data <- data[, c('ProteinName', 'FullPeptideName', 'Charge', 'aggr_Fragment_Annotation', 'aggr_Peak_Area', 
                     "RT","BioReplicate", "Condition", "Run")]

  colnames(data) <- gsub("Run", "filename", colnames(data))
  
  if(isTRUE(replace.Unimod)){
    # replace UniMod: to UniMod_
    data$FullPeptideName <- gsub("UniMod:", "UniMod_", data$FullPeptideName)
    data$aggr_Fragment_Annotation <- gsub("UniMod:", "UniMod_", data$aggr_Fragment_Annotation)
  }
    
  return(data)
}