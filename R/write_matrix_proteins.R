#' Writes out an overview matrix of summed signals per protein identifier
#' (lines) over run_id(columns).
#'
#' Writes out an overview matrix on protein level of a supplied (unfiltered or
#' filtered) OpenSWATH results data frame. The protein quantification is achieved
#' by summing the areas under all 6 transitions per precursor, summing all
#' precursors per FullPeptideName and all FullPeptideName signals per
#' ProteinName entry.
#' This function does not select consistently quantified or top peptides but
#' sums all signals availabe that may or may not originate from the same set of
#' peptides across different runs. A more detailed overview can be generated
#' using the function write_matrix_peptides().
#' Peptide selection can be achieved upstream using e.g. the functions
#' filter_mscore_requant(), filter_on_max_peptides() and
#' filter_on_min_peptides().
#'
#' @param data A data frame containing annotated OpenSWATH/pyProphet data.
#' @param write.csv  Option to determine if table should be written automatically into csv file.
#' @param fun.aggregate  What function to use when aggregating the set of intensities?
#' @param filename  File base name of the .csv matrix written out to the working
#'   folder
#' @param rm.decoy  Logical whether decoys will be removed from the data
#'   matrix. Defaults to FALSE. It's sometimes useful to know how decoys behave
#'   across a dataset and how many you allow into your final table with the
#'   current filtering strategy.
#' @return the peptides as a matrix, also output .csv matrix is written to the
#'   working folder
#' @author Moritz Heusel
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  written <- write_matrix_proteins(data)
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
