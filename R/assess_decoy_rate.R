#' Assess decoy rate in data
#'
#' This function assesses the number of quantifications (typically peptides)
#' that are decoys (false-positive) versus true identifications. 
#'
#' A printout is generated to indicate the number of non-decoy, decoy peptides
#' and the rate of decoy vs non-decoy peptides. Unique peptides are counted, so
#' a precursor with different charge states is counted as one peptide. In the
#' column "decoy" the values need to be 1,0 or TRUE and FALSE.
#'
#' @param data A data frame that contains at least a column named "FullPeptideName" and "decoy".
#' @param column.ids The column name of the Peptide identifier. Default: FullPeptideName.
#' @param column.decoy The column name of the decoy column. Default: decoy.
#' @return Message detailing the number of decoys, non-decoys, and the ratio.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data <- OpenSWATH_data
#'  assess_decoy_rate(data)
#' @export

assess_decoy_rate <- function(data, 
                              column.ids = "FullPeptideName",
                              column.decoy = "decoy") {
    if (sum(colnames(data) == column.decoy) < 1) {
        stop("There is no decoy column in the table")
    }
    if (sum(colnames(data) == column.ids) < 1) {
        stop("There is no column in the table for identifiers (e.g. peptides).")
    }

    add.colnames <- colnames(data)
    add.colnames <- add.colnames[add.colnames != column.decoy]

    .non_decoy.peptides <- unique(data[data$decoy == FALSE, c(column.ids)])
    .decoy.peptides <- unique(data[data$decoy == TRUE, c(column.ids)])

    message("Number of non-decoy identifiers: ", length(.non_decoy.peptides), "\n",
        "Number of decoy identifiers: ", length(.decoy.peptides), "\n", 
        "Decoy rate: ", sprintf("%.4f", (length(.decoy.peptides)/length(.non_decoy.peptides))))
}
