transform_MSstats_OpenSWATH <- function(data){
  colnames(data) <- gsub("\\.", "", colnames(data))

  missing_columns <- setdiff(c("ProteinName", 'PeptideSequence', 'PrecursorCharge', 'FragmentIon',
                                                         'ProductCharge', "IsotopeLabelType", 'Condition', "BioReplicate",
                                                         "Run", 'Intensity'), colnames(data))

  if(length(missing_columns > 0)){
    stop(paste("The data frame doesn't contain all required columns. The following columns are missing:\n",
               paste(missing_columns, collapse=", ")))

  }

  add.colnames <- colnames(data)[!(colnames(data) %in% c("ProteinName", 'PeptideSequence', 'PrecursorCharge', 'FragmentIon',
                                                         'ProductCharge', "IsotopeLabelType", 'Condition', "BioReplicate",
                                                         "Run", 'Intensity'))]

  data <- data[,c("ProteinName", 'PeptideSequence', 'PrecursorCharge', 'FragmentIon',
                  'ProductCharge', "IsotopeLabelType", 'Condition', "BioReplicate",
                  "Run", 'Intensity', add.colnames)]
  colnames(data)[seq_len(10)] <- c("ProteinName", 'FullPeptideName', 'Charge', 'FragmentIon',
                  'ProductCharge', "IsotopeLabelType", 'Condition', "BioReplicate",
                  "Run", 'Intensity')

  if(length(add.colnames) > 0){
    message("Additional columns present in the data: ", paste(add.colnames, collapse=", "))
  }

  if(!("mscore" %in% colnames(data))){
    message("No column 'mscore' present in the data. This column is required for functions estimating FDR or filtering the data.")
  }
  return(data)

}