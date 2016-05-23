sample_annotation <- function(data, sample.annotation, data.type="openSWATH", column.file = "align_origfilename", change.run.id = TRUE, verbose=FALSE){
  #### annotate sample
  ### needs a txt file with the columns Filename, Condition, BioReplicate, Run. In Filename a unique string contained in File
  ### must be contained.
  if(!(column.file %in% colnames(data))){
    warning("Warning: column for filename is not present in data file")
  }
  if(nlevels(factor(sample.annotation$Filename)) != nlevels(factor(data[,column.file]))){
    stop("The number of sample annotation condition and filenames in data are not balanced.", "\n",
         "Different filenames in sample annotation file: ", nlevels(factor(sample.annotation$Filename)), "\n",
         "Different filenames in data file: ", nlevels(factor(data[,column.file])))
  }
  
  for(i in 1:nrow(sample.annotation)){
    n.found <- grep(sample.annotation[i,"Filename"], sample.annotation[,"Filename"])
    if(length(n.found) > 1){
      stop("The values in the column filename are not unique and will lead to erroneous results because the following string matches to multiple different filenames: ", sample.annotation[i,"Filename"])
    }
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
                                                           'aggr_Peak_Area', 'Condition', "BioReplicate", "Run"))]

    data <- data[,c('ProteinName', 'FullPeptideName', 'Charge', 'aggr_Fragment_Annotation',
                    'aggr_Peak_Area', 'Condition', "BioReplicate", "Run", add.colnames)]
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

