filter_mscore <- function(data, mscore, rm.decoy=TRUE){
  if(sum(colnames(data) == "decoy") == 1 & rm.decoy == TRUE){
    data <- data[data$decoy == 0,]
    #subset(data, decoy == 0)
  }
  
  #data.filtered <- subset(data, m_score <= mscore)
  data.filtered <- data[data$m_score <= mscore,]
  
  message("Dimension difference: ", paste(dim(data)-dim(data.filtered), collapse=", "))
  
  return(data.filtered)
}

