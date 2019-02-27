filter_mscore <- function(data, mscore, rm.decoy=TRUE, mscore.col = "m_score"){
  if(sum(colnames(data) == "decoy") == 1 & rm.decoy == TRUE){
    data <- data[data$decoy == 0,]
    #subset(data, decoy == 0)
  }
  
  mscore.col <- JPP_update(data, mscore.col)

  data.filtered <- data[data[,mscore.col] <= mscore,]
  
  message("Dimension difference: ", paste(dim(data)-dim(data.filtered), collapse=", "))
  
  return(data.filtered)
}

