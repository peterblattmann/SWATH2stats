#' Query a swath2stats data set for decoy hit rates.
#'
#' This looks at the swath2stats table and compares the length of the detected
#' decoy peptides vs non-decoy peptides.
#'
#' @param data SWATH2stats data structure.
#' @return list containing the number of decoys, non-decoys, and the ratio.
#' @export
assess_decoy_rate <- function(data) {
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
