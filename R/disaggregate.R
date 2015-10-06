disaggregate <- function(data){
  
  data.new <- cbind(data, colsplit(data$aggr_Fragment_Annotation, ";", c("Split_FragAnnot_1", "Split_FragAnnot_2", "Split_FragAnnot_3", 
                                                             "Split_FragAnnot_4", "Split_FragAnnot_5", "Split_FragAnnot_6")),
                    colsplit(data$aggr_Peak_Area, ";", c("Split_PeakArea_1", "Split_PeakArea_2", "Split_PeakArea_3", 
                                                                   "Split_PeakArea_4", "Split_PeakArea_5", "Split_PeakArea_6")), stringsAsFactors=FALSE)
  
  data.new.m <- reshape2::melt(data.new, id.vars=grep("Split_FragAnnot", colnames(data.new), invert=TRUE), measure.vars=grep("Split_FragAnnot", colnames(data.new)), 
                     variable.name="FragAnnot_N", value.name = "Fragment")
  data.new.m2 <- reshape2::melt(data.new.m, id.vars=grep("Split_PeakArea", colnames(data.new.m), invert=TRUE), measure.vars=grep("Split_PeakArea", colnames(data.new.m)), 
                     variable.name="Area_N", value.name = "Area")
  
  # added because it didn't name the variables in a later trial
  if(sum(colnames(data.new.m2) %in% c("FragAnnot_N"))==0){
    l <- length(colnames(data.new.m2))
    colnames(data.new.m2)[l-3] <- "FragAnnot_N"
    colnames(data.new.m2)[l-2] <- "Fragment"
    colnames(data.new.m2)[l-1] <- "Area_N"
    colnames(data.new.m2)[l] <- "Area"
  }
  
  
  data.new.m3 <- data.new.m2[gsub("Split_FragAnnot_", "", data.new.m2[,"FragAnnot_N"]) == gsub("Split_PeakArea_", "", data.new.m2[,"Area_N"]),]
  
  cols <- colnames(data.new.m3)[colnames(data.new.m3) %in% c('ProteinName', 'FullPeptideName', "PeptideSequence", 'Sequence',
                                                             'Charge', "PrecursorCharge",'Fragment', "FragmentIon", 
                                                             'Area', 'Condition', "BioReplicate", "Run", "RT")]

  data.new.merged <- data.new.m3[,cols]
  
  colnames(data.new.merged) <- gsub("FullPeptideName", "PeptideSequence", colnames(data.new.merged))
  colnames(data.new.merged) <- gsub("^Charge$", "PrecursorCharge", colnames(data.new.merged))
  
  colnames(data.new.merged) <- gsub("Area", "Intensity", colnames(data.new.merged))
  colnames(data.new.merged) <- gsub("Fragment", "FragmentIon", colnames(data.new.merged))
  
  if("Sequence" %in% cols){
    colnames(data.new.merged) <- gsub("^Sequence$", "NakedSequence", colnames(data.new.merged))
  }
  
  
  return(data.new.merged)
}

