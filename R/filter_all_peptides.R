#' Select all proteins that are supported by peptides.
#'
#' This functions counts all proteins that are supported by peptides (including
#' non proteo-typic peptides). All peptides (incl. non proteotypic peptides are
#' selected. For the proteins supproted by proteotypic peptide the "1/" in front
#' of the identifier is removed to facilitate further data processing.
#'
#' @param data  A data frame containing SWATH data.
#' @param column  Which column contains the data to modify?
#' @param n  How many new IDs should we print to show if this worked?
#' @return  Returns a data frame with the data from both proteotypic and
#'   non-proteotypic peptides.
#' @author Peter Blattmann
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.all <- filter_all_peptides(data.filtered.decoy)
#' }
#' @export
filter_all_peptides <- function(data, column="proteinname", n=6) {
  data_proteins <- gsub(pattern="^[[:digit:]]\\/", replacement="", x=data[[column]])
  data_proteins <- unlist(strsplit(as.character(data_proteins), "\\/"))
  data_proteins <- unique(data_proteins)
  message("Number of proteins detected: ", length(data_proteins))
  message("First ", n, " protein identifiers: ", paste(head(data_proteins, n=n), collapse=", "))
  ## Delete non-unique features (unique features have "1/" in front of protein name)
  data[[column]] <- gsub(pattern="^1\\/", replacement="", x=data[[column]])
  return(data)
}
