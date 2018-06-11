write_matrix_proteins <- function(data, write.csv=FALSE,
                                  filename="SWATH2stats_overview_matrix_proteinlevel.csv",
                                  rm.decoy=FALSE) {
  if (rm.decoy == TRUE) {
    data <- subset(data, data[["decoy"]] == 0)
  }

  data.protein.sumall <- aggregate(data[,"intensity"], by=list(data[["proteinname"]], data[["run_id"]]), sum)
  colnames(data.protein.sumall) <- c("proteinname", "run_id", "intensity.all.sum")
  data.protein.sum.table <- reshape2::dcast(data.protein.sumall, proteinname ~ run_id,
                                            value.var="intensity.all.sum",
                                            fun.aggregate=sum)
  if (isTRUE(write.csv)) {
    write.csv(data.protein.sum.table, file=filename, row.names=FALSE, quote=FALSE)
    message("Protein overview matrix ", filename," written to working folder.")
  }
  return(data.protein.sum.table)
}
