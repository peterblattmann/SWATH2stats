filter_mscore <- function(data, mscore, rm.decoy=TRUE) {
  if (sum(colnames(data) == "decoy") == 1 & rm.decoy == TRUE) {
    data <- data[data[["decoy"]] == 0, ]
  }

  data.filtered <- data[data[["m_score"]] <= mscore, ]
  message("Original dimension: ", nrow(data), ", new dimension: ", nrow(data.filtered),
          ", difference: ", nrow(data) - nrow(data.filtered), ".")
  return(data.filtered)
}
