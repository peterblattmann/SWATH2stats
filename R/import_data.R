import_data <- function(data){
  colnames_OpenSWATH <- c('ProteinName', 'FullPeptideName', 'Charge', 'Sequence', 'aggr_Fragment_Annotation',
                          'aggr_Peak_Area', "transition_group_id", "decoy", "m_score", "RT", "align_origfilename", 
                          "Intensity", "not applicable")
    
  message("When reading in data, please specify to which column in the OpenSWATH data they correspond. 
For columns that do not correspond to any OpenSWATH column, choose \"not applicable\". 
Explanation of the OpenSWATH columns can be found in the manual page.\n"
)  
  
  # Dialogue to map columns
  for(i in colnames(data)){
    value <- select.list(colnames_OpenSWATH, title = paste("Select column that", i, "corresponds to the column."))
    if(value != "not applicable"){
      colnames(data) <- gsub(i, value, colnames(data))
      message("Column name ", i, " was substituted by ", value, "\n")
    }
    
  }
  
  # add NA to colnames that were not mapped
  add.colnames <- colnames(data)[!(colnames(data) %in% colnames_OpenSWATH)]
  if(length(add.colnames > 0)){
    for(i in add.colnames)
    data[,i] <- NA
    warning("Not all columns required within SWATH2stats were mapped: The columns ", paste(add.colnames, collapse =", "), " were added with NA as value.")
  }
  
  return(data)
}
