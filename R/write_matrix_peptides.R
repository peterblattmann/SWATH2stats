write_matrix_peptides<- function(data, write.csv=FALSE, filename = "SWATH2stats_overview_matrix_peptidelevel.csv", rm.decoy = FALSE)
 {
  if(rm.decoy == TRUE){
    data<-subset(data, data$decoy == 0)
  }
  data.peptide <- data[,c("ProteinName", "run_id", "FullPeptideName", "Intensity")]
  ProteinName_FullPeptideName<-paste(data.peptide$ProteinName,data.peptide$FullPeptideName, sep="_")
  data.peptide<-cbind(ProteinName_FullPeptideName, data.peptide)
  data.peptide.table<-dcast(data.peptide, ProteinName_FullPeptideName~run_id, value.var = "Intensity", fun.aggregate=sum)
  if(isTRUE(write.csv)){
    write.csv(data.peptide.table, file=filename, row.names=FALSE, quote=FALSE)
    message("Peptide overview matrix ", filename," written to working folder." , "\n")
  }
  return(data.peptide.table)
}