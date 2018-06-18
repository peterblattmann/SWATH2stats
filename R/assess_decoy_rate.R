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
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data <- OpenSWATH_data
#'  assess_decoy_rate(data)
#' }
#' @export
assess_decoy_rate <- function(data) {
  colnames(data) <- tolower(colnames(data))
  if(sum(colnames(data) == "decoy") < 1) {
    stop("There is no decoy column in the table")
  }
  if(sum(colnames(data) == "fullpeptidename") < 1) {
    stop("There is no fullpeptidename column in the table")
  }

  add_colnames <- colnames(data)
  add_colnames <- add_colnames[add_colnames != "decoy"]

  num_non_decoy_peptides <- unique(data[data[["decoy"]] == FALSE, c("fullpeptidename")])
  num_decoy_peptides <- unique(data[data[["decoy"]] == TRUE, c("fullpeptidename")])
  decoy_rate <- sprintf("%.4f", (length(num_decoy_peptides) / length(num_non_decoy_peptides)))

  message(
    "Number of non-decoy peptides: ", length(num_non_decoy_peptides), "\n",
    "Number of decoy peptides: ", length(num_decoy_peptides), "\n",
    "Decoy rate: ", decoy_rate
  )
  retlist <- list(
    "non_decoys" = num_non_decoy_peptides,
    "decoys" = num_decoy_peptides,
    "ratio" = decoy_rate)
  return(retlist)
}
