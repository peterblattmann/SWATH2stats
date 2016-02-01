sample_annotation <- function(data, sample.annotation, data.type="openSWATH", column.file = "align_origfilename", column.runid = "run_id", change.run.id = TRUE, verbose=FALSE){
  #### annotate sample
  ### needs a txt file with the columns Filename, Condition, BioReplicate, Run. In Filename a unique string contained in File
  ### must be contained.
  if(!(column.file %in% colnames(data))){
    warning("Warning: column for filename is not present in data file")
  }
  if(nlevels(factor(paste(sample.annotation$Filename))) != nlevels(factor(data[,column.file]))){
    stop("Warning: the number of sample annotation condition and filenames in data are not balanced.", "\n",
         "Different filenames in sample annotation file: ", nlevels(factor(sample.annotation$Filename)), "\n",
         "Different filenames in data file: ", nlevels(factor(data[,column.file])))
  }
  if(data.type=="openSWATH"){
    for(i in levels(factor(sample.annotation$Filename))){
      if(verbose){
        print(i)
      }
      coord <- grep(i, data[,column.file])

      if(length(coord) == 0){
        warning("No measurement value found for this sample in the data file: ", print(i))
      }
      data.subset <- sample.annotation[which(i == sample.annotation$Filename),]
      data[coord, "Condition"] <- data.subset[, "Condition"]
      data[coord, "BioReplicate"] <- data.subset[, "BioReplicate"]
      data[coord, "Run"] <- data.subset[, "Run"]
    }

    add.colnames <- colnames(data)[!(colnames(data) %in% c('ProteinName', 'FullPeptideName', 'Charge', 'aggr_Fragment_Annotation',
                                                           'aggr_Peak_Area', 'Condition', "BioReplicate", "Run", column.runid))]

    data <- data[,c('ProteinName', 'FullPeptideName', 'Charge', 'aggr_Fragment_Annotation',
                    'aggr_Peak_Area', 'Condition', "BioReplicate", "Run", column.runid, add.colnames)]
    if(change.run.id){
      data$run_id <- paste(data$Condition, data$BioReplicate, data$Run, sep="_")
    }
    return(data)
  }
  if(data.type=="MSstats"){
    colnames(data) <- gsub("Run", column.file, colnames(data))
    for(i in levels(sample.annotation$Filename)){
      if(verbose){
        print(i)
      }
      coord <- grep(i, data[,column.file])
      if(length(coord) == 0){
        warning("No measurement value found for this sample in the data file: ", print(i))
      }
      data.subset <- sample.annotation[which(i == sample.annotation$Filename),]
      data[coord, "Condition"] <- data.subset[, "Condition"]
      data[coord, "BioReplicate"] <- data.subset[, "BioReplicate"]
      data[coord, "Run"] <- data.subset[, "Run"]
    }

    add.colnames <- colnames(data)[!(colnames(data) %in% c("ProteinName", 'PeptideSequence', 'PrecursorCharge', 'FragmentIon',
                                                           'ProductCharge', "IsotopeLabelType", 'Condition', "BioReplicate",
                                                           "Run", 'Intensity'))]

    data <- data[,c("ProteinName", 'PeptideSequence', 'PrecursorCharge', 'FragmentIon',
                    'ProductCharge', "IsotopeLabelType", 'Condition', "BioReplicate",
                    "Run", 'Intensity', add.colnames)]

    return(data)
  }
  else(
    print("Type of data not recognized")
  )
}

