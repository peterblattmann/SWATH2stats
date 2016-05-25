utils::globalVariables(c("aggregate", "head"))

convert4mapDIA <- function(data, RT=FALSE){

  mDIA.columns <- c('ProteinName', 'PeptideSequence', 'FragmentIon',
                       "Intensity", "BioReplicate", "Condition", "RT")

    col.names.missing <- mDIA.columns[!(mDIA.columns %in% colnames(data))]
  if(length(col.names.missing) > 0){
    warning("One or several columns required by mDIA were not in the data and filled with NAs.
            Missing columns: ", paste(unlist(col.names.missing), collapse=", "))

    data[, col.names.missing] <- NA
  }

  # replace UniMod: to UniMod_
  data$PeptideSequence <- gsub(":", "_", data$PeptideSequence)
  data$FragmentIon     <- gsub(":", "_", data$FragmentIon)

  data.red.wide.test <- dcast(data, ProteinName + PeptideSequence + FragmentIon ~ Condition + BioReplicate, fun.aggregate=length, value.var="Intensity")

  if(sum(data.red.wide.test[,4:dim(data.red.wide.test)[2]] > 1)>0){
    col.names.repl <- apply(data.red.wide.test[,4:dim(data.red.wide.test)[2]], 2, function(x) sum(x > 1))
    col.names.repl <- names(col.names.repl[col.names.repl > 0])
    row.names.repl <- data.red.wide.test[apply(data.red.wide.test[,4:dim(data.red.wide.test)[2]], 1, function(x) sum(x > 1)), "FragmentIon"]
    
    warning("Data contains several intensity values per condition\n\n",
            "in the following columns: ", paste(col.names.repl, collapse=", "),"\n\n",
            "and in the following rows: ", paste(row.names.repl, collapse=", "))
  }

  data.red.wide <- dcast(data, ProteinName + PeptideSequence + FragmentIon ~ Condition + BioReplicate, fun.aggregate=mean, value.var="Intensity")
  if(RT){
    RTs <- unique(aggregate(data[,c(which(colnames(data)=="RT"))],
                            by=list(data$FragmentIon), FUN=mean, na.rm=TRUE))
    colnames(RTs) <- c("FragmentIon", "RT")
    RTs$RT <- RTs$RT / 60
    data.red.wide <- merge(data.red.wide, RTs, by="FragmentIon", all.x=TRUE)
    data.red.wide <- data.red.wide[,c("ProteinName", "PeptideSequence", "FragmentIon",
                                      colnames(data.red.wide)[!(colnames(data.red.wide) %in% c("ProteinName", "PeptideSequence", "FragmentIon", "RT"))],
                                      "RT")]

  }

  return(data.red.wide)
}
