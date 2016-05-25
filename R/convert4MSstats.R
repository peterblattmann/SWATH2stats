convert4MSstats <- function(data, replace.values = TRUE, replace.colnames = TRUE, replace.Unimod = TRUE){

  MSstats.columns <- c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon',
                       'ProductCharge', "IsotopeLabelType", "Intensity", "BioReplicate", "Condition", "Run")

  if(isTRUE(replace.colnames)){
    colnames(data) <- gsub("FullPeptideName", "PeptideSequence", colnames(data))
    colnames(data) <- gsub("^Charge$", "PrecursorCharge", colnames(data))
    colnames(data) <- gsub("align_origfilename", "File", colnames(data))
    }
    col.names.missing <- MSstats.columns[!(MSstats.columns %in% colnames(data))]
    if(length(col.names.missing) > 0){
      message("One or several columns required by MSstats were not in the data. The columns were created and filled with NAs.\nMissing columns: ",
              paste(unlist(col.names.missing), collapse=", "))

      if("IsotopeLabelType" %in% col.names.missing){
        message("IsotopeLabelType was filled with light.")
      }

      data[, col.names.missing] <- NA

      if("IsotopeLabelType" %in% col.names.missing){
        data[, "IsotopeLabelType"] <- "light"
      }
      if("PrecursorCharge" %in% col.names.missing){
          data[, "PrecursorCharge"] <- gsub(".*_([[:digit:]])$", "\\1", data[, "FragmentIon"])
      }
    }
    data <- data[, MSstats.columns]

    if(isTRUE(replace.values)){
    # replace negative values to 0 and 0 to NA
    if(sum(data$Intensity < 0, na.rm=TRUE) > 0){
      data[data$Intensity < 0,"Intensity"] <- 0
      warning("Negative intensity values were replaced by NA")
    }
    if(sum(data$Intensity == 0, na.rm=TRUE) > 0){
      data[data$Intensity %in% 0,"Intensity"] <- NA
      warning("Intensity values that were 0, were replaced by NA")
    }
    }

    if(isTRUE(replace.Unimod)){
      # replace UniMod: to UniMod_
      data$PeptideSequence <- gsub("UniMod:", "UniMod_", data$PeptideSequence)
      data$FragmentIon     <- gsub("UniMod:", "UniMod_", data$FragmentIon)
    }
    return(data)
}
