convert4MSstats <- function(data, replace.values=TRUE,
                            replace.colnames=TRUE, replace.Unimod=TRUE) {

  MSstats.columns <- c("proteinname", "peptidesequence", "precursorcharge",
                       "fragmention", "productcharge", "isotopelabeltype",
                       "intensity", "bioreplicate", "condition", "run")

  if (isTRUE(replace.colnames)) {
    colnames(data) <- gsub("fullpeptidename", "peptidesequence", colnames(data))
    colnames(data) <- gsub("^charge$", "precursorcharge", colnames(data))
    colnames(data) <- gsub("align_origfilename", "file", colnames(data))
  }
  col.names.missing <- MSstats.columns[!(MSstats.columns %in% colnames(data))]

  if (length(col.names.missing) > 0) {
      message("One or several columns required by MSstats were not in the data. The columns were created and filled with NAs.\nMissing columns: ",
              paste(unlist(col.names.missing), collapse=", "))

      data[, col.names.missing] <- NA

      if ("isotopelabeltype" %in% col.names.missing) {
        message("isotopelabeltype was filled with light.")
        data[, "isotopelabeltype"] <- "light"
      }

      if ("precursorcharge" %in% col.names.missing) {
          data[, "precursorcharge"] <- gsub(".*_([[:digit:]])$", "\\1", data[, "fragmention"])
      }
    }
    data <- data[, MSstats.columns]

    if (isTRUE(replace.values)) {
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

    if(isTRUE(replace.Unimod)){
      ## replace UniMod: to UniMod_
      data[["peptidesequence"]] <- gsub("UniMod:", "UniMod_", data[["peptidesequence"]])
      data[["fragmention"]] <- gsub("UniMod:", "UniMod_", data[["fragmention"]])
    }
  return(data)
}
