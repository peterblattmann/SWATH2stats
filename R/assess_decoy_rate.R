assess_decoy_rate <- function(data) {
  if(sum(colnames(data) == "decoy") < 1) {
    stop("There is no decoy column in the table")
  }
  if(sum(colnames(data) == "fullpeptidename") < 1) {
    stop("There is no fullpeptidename column in the table")
  }

  add.colnames <- colnames(data)
  add.colnames <- add.colnames[add.colnames != "decoy"]

  .non_decoy.peptides <- unique(data[data[["decoy"]] == FALSE, c("fullpeptidename")])
  .decoy.peptides <- unique(data[data[["decoy"]] == TRUE, c("fullpeptidename")])

  message(
    "Number of non-decoy peptides: ", length(.non_decoy.peptides), "\n",
    "Number of decoy peptides: ", length(.decoy.peptides), "\n",
    "Decoy rate: ", sprintf("%.4f", (length(.decoy.peptides)/length(.non_decoy.peptides)))
  )
}
