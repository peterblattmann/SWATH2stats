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
#' @param column Column to query for decoy string
#' @param decoy_string String defining a decoy. Default: DECOY
#' @return Returns a data frame with changed protein labels
#' @author Moritz Heusel
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.2 <- removeDecoyProteins(data.filtered.decoy)
#' @export
removeDecoyProteins <- function(data,
                                column = "ProteinName",
                                decoy_string = "DECOY") {
    ids <- grep(decoy_string, data[,column])
    identifiers <- data[ids, column]
    identifiers_split <- strsplit(as.character(identifiers), "/")
    identifiers_split_removed <- lapply(identifiers_split, rmDecoyProt)
    identifiers_removed <- vapply(identifiers_split_removed, function(x) {
        paste(x, collapse = "/")
    }, "a")
    data[ids, "ProteinName"] <- identifiers_removed
    return(data)
}

#' Subfunction to remove decoys
#'
#' @param x proteinname string to query.
#' @param decoy_string String defining a decoy
#' @return returns string without elements containing the decoy string
rmDecoyProt <- function(x,
                        decoy_string = "DECOY") {
    ids <- grep(decoy_string, x)
    ids.sel <- grep(decoy_string, x, invert = TRUE)
    x[1] <- as.numeric(x[1]) - length(ids)
    x <- x[ids.sel]
    return(x)
}
