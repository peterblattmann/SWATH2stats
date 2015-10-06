write_matrix_proteins<- function(data, filename = "SWATH2stats_overview_matrix_proteinlevel.csv", rm.decoy = FALSE)
{
  if(rm.decoy == TRUE){
    data<-subset(data, data$decoy == 0)
  }
  data.protein.sumall <- aggregate(data[,"Intensity"], by=list(data$ProteinName, data$run_id), sum)
  colnames(data.protein.sumall) <- c("ProteinName", "run_id", "Intensity.all.sum")
  data.protein.sum.table <- dcast(data.protein.sumall, ProteinName~run_id, fun.aggregate=sum)
  write.csv(data.protein.sum.table, file=filename, row.names=FALSE, quote=FALSE)
  message("Protein overview matrix ", filename," written to working folder.", "\n" )
}
