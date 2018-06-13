#' Take data from the swath2stats format and get it ready for use by MSstats.
#'
#' Though the two tools use very similar formats, some coercion is required to
#' pass swath2stats data into msstats; ergo this function.
#'
#' @param data  SWATH2stats formatted data.
#' @param replace_values  Change the peculiar sample IDs from openswath into the
#'   somewhat more human-readable format generated from a combination of
#'   run/replicate/etc.
#' @param replace_colnames  Change a few column names to make them fit better in
#'   msstats?
#' @param replace_unimod  Change the unimod column name to make it fit better in
#'   msstats?
#' @return MSstats compatible data, hopefully.
#' @export
convert_for_msstats <- function(data, replace_values=TRUE,
                                replace_colnames=TRUE, replace_unimod=TRUE) {

  msstats_columns <- c("proteinname", "peptidesequence", "precursorcharge",
                       "fragmention", "productcharge", "isotopelabeltype",
                       "intensity", "bioreplicate", "condition", "run")

  if (isTRUE(replace_colnames)) {
    colnames(data) <- gsub(pattern="fullpeptidename",
                           replacement="peptidesequence",
                           x=colnames(data))
    colnames(data) <- gsub(pattern="^charge$", replacement="precursorcharge",
                           x=colnames(data))
    colnames(data) <- gsub(pattern="align_origfilename", replacement="file",
                           x=colnames(data))
  }
  col_names_missing <- msstats_columns[!(msstats_columns %in% colnames(data))]

  if (length(col_names_missing) > 0) {
      message("One or several columns required by MSstats were not in the data. The columns were created and filled with NAs.\nMissing columns: ",
              paste(unlist(col_names_missing), collapse=", "))

      data[, col_names_missing] <- NA

      if ("isotopelabeltype" %in% col_names_missing) {
        message("isotopelabeltype was filled with light.")
        data[, "isotopelabeltype"] <- "light"
      }

      if ("precursorcharge" %in% col_names_missing) {
          data[, "precursorcharge"] <- gsub(".*_([[:digit:]])$", "\\1", data[, "fragmention"])
      }
    }
    data <- data[, msstats_columns]

    if (isTRUE(replace_values)) {
    # replace negative values to 0 and 0 to NA
      if (sum(data[["intensity"]] < 0, na.rm=TRUE) > 0) {
        data[data[["intensity"]] < 0, "intensity"] <- 0
        warning("Negative intensity values were replaced by NA")
      }

      if (sum(data[["intensity"]] == 0, na.rm=TRUE) > 0) {
        data[data[["intensity"]] %in% 0, "intensity"] <- NA
      warning("Intensity values that were 0, were replaced by NA")
      }
    }

    if(isTRUE(replace_unimod)){
      ## replace UniMod: to UniMod_
      data[["peptidesequence"]] <- gsub("UniMod:", "UniMod_", data[["peptidesequence"]])
      data[["fragmention"]] <- gsub("UniMod:", "UniMod_", data[["fragmention"]])
    }
  return(data)
}
