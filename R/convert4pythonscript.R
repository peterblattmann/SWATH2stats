convert4pythonscript <- function(data, replace.Unimod=TRUE) {

  data <- data[, c("proteinname", "fullpeptidename", "charge", "aggr_fragment_annotation", "aggr_peak_area",
                   "rt","bioreplicate", "condition", "run")]

  colnames(data) <- gsub("run", "filename", colnames(data))

  if (isTRUE(replace.Unimod)) {
    # replace UniMod: to UniMod_
    data[["fullpeptidename"]] <- gsub("UniMod:", "UniMod_", data[["fullpeptidename"]])
    data[["aggr_fragment_annotation"]] <- gsub("UniMod:", "UniMod_", data[["aggr_fragment_annotation"]])
  }

  return(data)
}
