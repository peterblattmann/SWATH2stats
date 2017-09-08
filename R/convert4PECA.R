utils::globalVariables(c("aggregate", "head"))

convert4PECA <- function(data){

  columns <- c('ProteinName', 'FullPeptideName', "Charge", "Intensity", "BioReplicate", "Condition", "RT")

    col.names.missing <- columns[!(columns %in% colnames(data))]
  if(length(col.names.missing) > 0){
    warning("One or several columns required by mDIA were not in the data and filled with NAs.
            Missing columns: ", paste(unlist(col.names.missing), collapse=", "))

    data[, col.names.missing] <- NA
  }

  data$FullPeptideName_Charge <- paste(data$FullPeptideName, data$Charge, sep= "_")
    
  # # replace UniMod: to UniMod_
  # data$PeptideSequence <- gsub(":", "_", data$PeptideSequence)
  # data$FragmentIon     <- gsub(":", "_", data$FragmentIon)
  
  data.red.wide.test <- dcast(data, ProteinName + FullPeptideName_Charge ~ Condition + BioReplicate, fun.aggregate=length, value.var="Intensity")
  
  if(sum(data.red.wide.test[,4:dim(data.red.wide.test)[2]] > 1)>0){
    col.names.repl <- apply(data.red.wide.test[,4:dim(data.red.wide.test)[2]], 2, function(x) sum(x > 1))
    col.names.repl <- names(col.names.repl[col.names.repl > 0])
    row.names.repl <- apply(data.red.wide.test[,4:dim(data.red.wide.test)[2]], 1, function(x) sum(x > 1))
    row.names.repl <- data.red.wide.test[row.names.repl > 0, "FullPeptideName_Charge"]
    
    warning("Data contains several intensity values per condition\n\n",
            "in the following columns: ", paste(col.names.repl, collapse=", "),"\n\n",
            "and in the following rows: ", paste(row.names.repl, collapse=", "))
  }
  
  data.out <- dcast(data, ProteinName + FullPeptideName_Charge  ~ Condition + BioReplicate, fun.aggregate=mean, value.var="Intensity")
  
  data.out <- droplevels(data.out)
  
  return(data.out)
}
