#' Removes decoy proteins from the protein group label
#'
#' There exist peptides annotated as protein groups with
#' 2/ProteinA/DECOY_ProteinB. However these are in principal proteotypic
#' peptides and should be annoated 1/ProteinA. This function changes these
#' labels accordingly. The subfunction rmDecoyProt removes the Decoy protein,
#' calling removeDecoyProteins also changes the nubmer before the protein group
#' accordingly.
#'
#' @param data A data frame containing SWATH data.
#' @param column which column to query?
#' @return Returns a data frame with changed protein labels
#' @author Moritz Heusel
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.2 <- removeDecoyProteins(data.filtered.decoy)
#' }
#' @export
removeDecoyProteins <- function(data, column="proteinname") {
  ids <- grep("^[1-9][0-9]*[0-9]*/.*DECOY.*", data[[column]])
  identifiers <- data[ids, column]
  identifiers_split <- strsplit(as.character(identifiers), "/")
  identifiers_split_removed <- lapply(identifiers_split, rmDecoyProt)
  identifiers_removed <- vapply(identifiers_split_removed, function(x) { paste(x, collapse="/")}, "a")
  data[ids, column] <- identifiers_removed
  return(data)
}

#' What string defines decoys?
#'
#' @param x proteinname string to query.
#' @param pattern chosen string to seek out.
rmDecoyProt <- function(x, pattern="DECOY") {
  ids <- grep(pattern=pattern, x=x)
  ids.sel <- grep(pattern=pattern, x=x, invert=TRUE)
  x[1] <- as.numeric(x[1]) - length(ids)
  x <- x[ids.sel]
  return(x)
}
