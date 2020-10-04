#' Select all proteins that are supported by peptides.
#'
#' This function can be used as opposed to the function 
#' "filter_proteotypic_peptides()". This function counts all proteins (including
#' proteins supported by non proteo-typic (i.e. shared) peptides). 
#' All peptides (incl. non proteotypic peptides are selected. 
#' For the proteins supported by proteotypic peptide the "1/" in front
#' of the identifier is removed to facilitate further data processing. 
#' The protein identifier of shared peptides needs to be separated by a slash "/".
#'
#' @param data  A data frame containing SWATH data.
#' @return  Returns a data frame with the data from both proteotypic and
#'   non-proteotypic peptides.
#' @author Peter Blattmann
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.all <- filter_all_peptides(data.filtered.decoy)
#'  }
#' @importFrom utils head
#' @export
filter_all_peptides <- function(data) {

    data.proteins <- gsub("^[[:digit:]]\\/", "", data$ProteinName)
    data.proteins <- unlist(strsplit(as.character(data.proteins), "\\/"))
    data.proteins <- unique(data.proteins)

    message("Number of proteins detected: ", length(data.proteins), 
            "\n", "Protein identifiers: ",
        paste(head(data.proteins), collapse = ", "))

    # Delete non-unique features (unique features have '1/' in front of protein name)
    data$ProteinName <- gsub("^1\\/", "", data$ProteinName)

    return(data)
}
