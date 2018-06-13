#' Write out a matrix of proteins from swath2stats data.
#'
#' This function should be folded into the write_matrix_peptides.
#'
#' @param data SWATH2stats data to write out.
#' @param write.csv  Write the result to a csv file?  This should be folded into the next arg.
#' @param fun.aggregate  What function to use when aggregating the set of intensities?
#' @param filename  Where to write the data?
#' @param rm.decoy  Remove the decoys?
#' @return the peptides as a matrix!
#' @export
write_matrix_proteins <- function(data, write.csv=FALSE, fun.aggregate=sum,
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
