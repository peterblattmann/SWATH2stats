#' Filter data to keep rows which make sense with respect to proteolytic ends.
#'
#' @param data  SWATH2stats data to filter.
#' @param rm.decoy  Drop decoys?
#' @param column  Which column to query for filtering?
#' @return  Filtered data!
#' @export
filter_proteotypic_peptides <- function(data, rm.decoy=TRUE, column="proteinname") {
  data <- unifyProteinGroupLabels(data, column=column)
  if (isTRUE(rm.decoy)) {
    data <- removeDecoyProteins(data)
  }

  data.proteins <- gsub("^[[:digit:]]\\/", "", data[[column]])
  data.proteins <- unlist(strsplit(as.character(data.proteins), "\\/"))
  data.proteins <- unique(data.proteins)

  message("Number of proteins detected: ", length(data.proteins), "\n",
          "Protein identifiers: ", paste(head(data.proteins), collapse=", "))

  # Delete non-unique features (unique features have "1/" in front of protein name)
  data.unique <- data[grep("^1/", data[[column]]), ]
  data.unique[["proteinname"]] <- gsub("^1\\/", "", data.unique[[column]])
  message("Number of proteins detected that are supported by a proteotypic peptide: ",
          length(unique(data.unique[[column]])), "\n",
          "Number of proteotypic peptides detected: ",
          length(unique(data.unique[["fullpeptidename"]])))
  return(data.unique)
}
