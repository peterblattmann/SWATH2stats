#' Import some data into swath2stats, presumably from openMS.
#'
#' @param data  I am guessing tric-formatted tsv data?
#' @return swath2stats data
#' @export
import_data <- function(data) {
  colnames_OpenSWATH <- c("proteinname", "fullpeptidename", "charge", "sequence",
                          "aggr_fragment_annotation", "aggr_peak_area",
                          "transition_group_id", "decoy", "m_score",
                          "rt", "align_origfilename", "intensity", "not applicable")

  colnames(data) <- tolower(colnames(data))
  message("When reading in data, please specify to which column in the OpenSWATH data they correspond.
For columns that do not correspond to any OpenSWATH column, choose \"not applicable\".
Explanation of the OpenSWATH columns can be found in the manual page.")

  # Dialogue to map columns
  for (i in colnames(data)) {
    value <- select.list(
      colnames_OpenSWATH,
      title=paste0("Select a column from the OpenSWATH output that the column ", i, "corresponds to."))
    if (value != "not applicable") {
      colnames(data) <- gsub(pattern=i, replacement=value, x=colnames(data))
      message("Column name ", i, " was replaced with ", value, ".")
    }
  }

  # add NA to colnames that were not mapped
  add.colnames <- colnames(data)[!(colnames(data) %in% colnames_OpenSWATH)]
  if (length(add.colnames > 0)) {
    for(i in add.colnames) {
      data[, i] <- NA
    }
    warning("Not all columns required within SWATH2stats were mapped: The columns ",
            paste(add.colnames, collapse =", "), " were added with NA as value.")
  }

  return(data)
}
