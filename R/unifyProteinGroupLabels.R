#' Unify the protein group labels.
#'
#' Unify the protein group labels (2/ProteinA/ProteinB and 2/ProteinB/ProteinA)
#' to one common label (e.g. 2/ProteinA/ProteinB)
#'
#' @param data  A data frame containing SWATH data.
#' @param column Which column to use for unifying the groups.
#' @return Returns a data frame with the unififed protein labels.
#' @author Moritz Heusel
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.unified <- unifyProteinGroupLabels(data.filtered.decoy)
#' }
#' @export
unifyProteinGroupLabels <- function(data, column="proteinname") {
  data[[column]] <- as.character(data[[column]])
  ids <- grep("^([2-9])|([1-9][0-9][0-9]*)/", data[[column]])
  identifiers <- data[ids, column]
  identifiers_split <- strsplit(as.character(identifiers), "/")
  identifiers_split_sorted <- lapply(identifiers_split, function(x) { sort(x) })
  identifiers_sorted <- vapply(identifiers_split_sorted, function(x) { paste(x, collapse="/") }, "a")
  data[ids, column] <- identifiers_sorted
  return(data)
}
