#' Filter for proteins that are supported by proteotypic peptides.
#'
#' Peptides can match to several proteins. With this function proteotypic
#' peptides, peptides that are only contained in one protein are
#' selected. Additionally the number of proteins are counted and printed.
#'
#' @param data  A data frame containing SWATH data.
#' @param rm.decoy  Option to remove the decoys during filtering.
#' @param column  Which column to query for filtering?
#' @return  Returns a data frame with only the data supported by proteotypic
#'   peptides.
#' @author  Peter Blattmann
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.all <- filter_proteotypic_peptides(data.filtered.decoy)
#' }
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
