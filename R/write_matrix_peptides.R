#' Writes out an overview matrix of peptides mapping to a FDR quality controlled
#' protein master list at controlled global peptide FDR quality.
#'
#' Writes out an overview matrix on peptide level of a supplied (unfiltered or
#' prefiltered) OpenSWATH results data frame.
#' The peptide quantification is achieved by summing the areas under all 6
#' transitions per precursor and summing all precursors per FullPeptideName.
#' In order to keep the peptide-to-protein association, the FullPeptideName is
#' joined with the ProteinName.
#'
#' @param data A data frame containing annotated OpenSWATH/pyProphet data.
#' @param write_csv Option to determine if table should be written automatically
#'   into csv file.
#' @param fun_aggregate  What function to use when aggregating the set of 
#'   intensities (sum or mean)?. Default: sum.
#' @param filename File base name of the .csv matrix written out to the working
#'   folder.
#' @param rm_decoy Logical whether decoys will be removed from the data
#'   matrix. Defaults to FALSE. It's sometimes useful to know how decoys behave
#'   across a dataset and how many you allow into your final table with the
#'   current filtering strategy.
#' @return the peptides as a matrix! also output .csv matrix is written to the
#'   working folder.
#' @author Moritz Heusel
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  written <- write_matrix_peptides(data)
#'  }
#' @importFrom reshape2 dcast
#' @export
write_matrix_peptides <- function(data, 
                                  write_csv = FALSE, 
                                  fun_aggregate = "sum",
                                  filename = "SWATH2stats_overview_matrix_peptidelevel.csv",
                                  rm_decoy = FALSE) {
    if (rm_decoy == TRUE) {
        data <- subset(data, data$decoy == 0)
    }
    data.peptide <- data[, c("ProteinName", "run_id", "FullPeptideName", "Intensity")]
    ProteinName_FullPeptideName <- paste(data.peptide$ProteinName, 
                                         data.peptide$FullPeptideName,
                                         sep = "_")
    data.peptide <- cbind(ProteinName_FullPeptideName, data.peptide)
    if(fun_aggregate == "sum"){
        data.peptide.table <- dcast(data.peptide, 
                                ProteinName_FullPeptideName ~ run_id,
                                value.var = "Intensity", fun.aggregate = sum)
    }
    if(fun_aggregate == "mean"){
        data.peptide.table <- dcast(data.peptide, 
                                    ProteinName_FullPeptideName ~ run_id,
                                    value.var = "Intensity", fun.aggregate = mean)
    }
    if(write_csv) {
        write.csv(data.peptide.table, file = filename, row.names = FALSE, quote = FALSE)
        message("Peptide overview matrix ", filename, " written to working folder.",
            "\n")
    }
    return(data.peptide.table)
}
