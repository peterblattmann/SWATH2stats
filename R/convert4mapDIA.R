#' Convert table into the format for mapDIA
#'
#' This functions selects the columns necessary for mapDIA.
#'
#' @param data  A data frame containing SWATH data.
#' @param RT   Option to export the retention times.
#' @return Returns a data frame in the appropriate format for mapDIA.
#' @note The table must not contain any technical replica, the intensity of
#'   technical replica is averaged. This function requires the package
#'   reshape2.
#' @references Teo, G., et al. (2015). "mapDIA: Preprocessing and statistical
#'   analysis of quantitative proteomics data from data independent acquisition
#'   mass spectrometry." J Proteomics 129: 108-120.
#' @author Peter Blattmann
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#'  data.mapDIA <- convert4mapDIA(raw, RT=TRUE)
#' }
#' @export
convert4mapDIA <- function(data, RT=FALSE) {

  mDIA.columns <- c("proteinname", "peptidesequence", "fragmention",
                    "intensity", "bioreplicate", "condition", "rt")

  col.names.missing <- mDIA.columns[!(mDIA.columns %in% colnames(data))]
  if (length(col.names.missing) > 0) {
    warning("One or several columns required by mDIA were not in the data and filled with NAs.
            Missing columns: ", paste(unlist(col.names.missing), collapse=", "))

    data[, col.names.missing] <- NA
  }

  # replace UniMod: to UniMod_
  data[["peptidesequence"]] <- gsub(":", "_", data[["peptidesequence"]])
  data[["fragmention"]] <- gsub(":", "_", data[["fragmention"]])

  data.red.wide.test <- dcast(
    data,
    proteinname + peptidesequence + fragmention ~ condition + bioreplicate,
    fun.aggregate=length,
    value.var="intensity")

  if (sum(data.red.wide.test[, 4:dim(data.red.wide.test)[2]] > 1) > 0) {
    col.names.repl <- apply(data.red.wide.test[, 4:dim(data.red.wide.test)[2]], 2, function(x) sum(x > 1))
    col.names.repl <- names(col.names.repl[col.names.repl > 0])
    row.names.repl <- apply(data.red.wide.test[, 4:dim(data.red.wide.test)[2]], 1, function(x) sum(x > 1))
    row.names.repl <- data.red.wide.test[row.names.repl > 0, "fragmention"]

    warning("Data contains several intensity values per condition\n\n",
            "in the following columns: ", paste(col.names.repl, collapse=", "), "\n\n",
            "and in the following rows: ", paste(row.names.repl, collapse=", "))
  }

  data.red.wide <- dcast(
    data,
    proteinname + peptidesequence + fragmention ~ condition + bioreplicate,
    fun.aggregate=mean,
    value.var="intensity")

  if (RT) {
    RTs <- unique(aggregate(data[, c(which(colnames(data)=="rt"))],
                            by=list(data[["fragmention"]]), FUN=mean, na.rm=TRUE))
    colnames(RTs) <- c("fragmention", "rt")
    RTs[["rt"]] <- RTs[["rt"]] / 60
    data.red.wide <- merge(data.red.wide, RTs, by="fragmention", all.x=TRUE)
    data.red.wide <- data.red.wide[, c("proteinname", "peptidesequence", "fragmention",
                                       colnames(data.red.wide)[!(colnames(data.red.wide) %in% c("proteinname", "peptidesequence", "fragmention", "rt"))],
                                       "rt")]
  }
  return(data.red.wide)
}
