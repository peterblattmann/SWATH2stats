#' Query a swath2stats data set for decoy hit rates.
#'
#' This looks at the swath2stats table and compares the length of the detected
#' decoy peptides vs non-decoy peptides.
#'
#' A printout is generated to indicate the number of non-decoy, decoy peptides
#' and the rate of decoy vs non-decoy peptides. Unique peptides are counted, so
#' a precursor with different charge states is counted as one peptide. In the
#' column "decoy" the values need to be 1,0 or TRUE and FALSE.
#'
#' @param data A data frame that contains at least a column named "FullPeptideName" and "decoy".
#' @return list containing the number of decoys, non-decoys, and the ratio.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data <- OpenSWATH_data
#'  assessed <- assess_decoy_rate(data)
#'  summary(assessed)
#' @export
assess_decoy_rate <- function(data) {
    if (sum(colnames(data) == "decoy") < 1) {
        stop("There is no decoy column in the table")
    }
    if (sum(colnames(data) == "FullPeptideName") < 1) {
        stop("There is no FullPeptideName column in the table")
    }

    add.colnames <- colnames(data)
    add.colnames <- add.colnames[add.colnames != "decoy"]

    .non_decoy.peptides <- unique(data[data$decoy == FALSE, c("FullPeptideName")])
    .decoy.peptides <- unique(data[data$decoy == TRUE, c("FullPeptideName")])

    message("Number of non-decoy peptides: ", length(.non_decoy.peptides), "\n",
        "Number of decoy peptides: ", length(.decoy.peptides), "\n", "Decoy rate: ",
        sprintf("%.4f", (length(.decoy.peptides)/length(.non_decoy.peptides))))
}
