write_matrix_peptides<- function(data, filename = "SWATH2stats_overview_matrix_peptidelevel.csv", rm.decoy = FALSE) 
 {
  if(rm.decoy == TRUE){
    data<-subset(data, data$decoy == 0)
  }
  data.peptide <- data[,c("ProteinName", "run_id", "FullPeptideName", "Intensity")]
  ProteinName_FullPeptideName<-paste(data.peptide$ProteinName,"_",data.peptide$FullPeptideName, sep="")
  data.peptide<-cbind(ProteinName_FullPeptideName, data.peptide)
  data.peptide.table<-dcast(data.peptide, ProteinName_FullPeptideName~run_id, fun.aggregate=sum)
  write.csv(data.peptide.table, file=filename, row.names=FALSE, quote=FALSE)
  message("Peptide overview matrix ", filename," written to working folder." , "\n")
}