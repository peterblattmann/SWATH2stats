#'  Convert table into the format for ROPECA
#'
#' This functions selects the columns necessary for ROPECA.
#'
#' @param data   A data frame containing SWATH data.
#' @return Returns a data frame in the appropriate format for ROPECA.
#' @note The table must not contain any technical replica, the intensity of
#'   technical replica is averaged. This function requires the package reshape2.
#' @author Peter Blattmann
#' @references Suomi, T. and Elo L.L. (2017). "Enhanced differential expression
#'   statistics for data-independent acquisition proteomics" Scientific Reports
#'   7, Article number: 5869.doi:10.1038/s41598-017-05949-y
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.PECA <- convert_PECA(data.filtered.decoy)
#' }
#' @export
convert_PECA <- function(data) {
  columns <- c("proteinname", "fullpeptidename", "charge", "intensity",
               "bioreplicate", "condition", "rt")

  col.names.missing <- columns[!(columns %in% colnames(data))]
  if (length(col.names.missing) > 0) {
    warning("One or several columns required by mDIA were not in the data and filled with NAs.
            Missing columns: ", paste(unlist(col.names.missing), collapse=", "))
    data[, col.names.missing] <- NA
  }

  data[["fullpeptidename_charge"]] <- paste(data[["fullpeptidename"]], data[["charge"]], sep= "_")

  ## # replace UniMod: to UniMod_
  ## data$PeptideSequence <- gsub(":", "_", data$PeptideSequence)
  ## data$FragmentIon     <- gsub(":", "_", data$FragmentIon)

  ## Interesting, this uses length() for aggregate, why is that I wonder?
  data.red.wide.test <- reshape2::dcast(data, proteinname + fullpeptidename_charge ~ condition + bioreplicate,
                                        fun.aggregate=length, value.var="intensity")

  if (sum(data.red.wide.test[, 4:dim(data.red.wide.test)[2]] > 1) > 0) {
    col.names.repl <- apply(data.red.wide.test[, 4:dim(data.red.wide.test)[2]], 2, function(x) sum(x > 1))
    col.names.repl <- names(col.names.repl[col.names.repl > 0])
    row.names.repl <- apply(data.red.wide.test[, 4:dim(data.red.wide.test)[2]], 1, function(x) sum(x > 1))
    row.names.repl <- data.red.wide.test[row.names.repl > 0, "fullpeptidename_charge"]

    warning("Data contains several intensity values per condition\n\n",
            "in the following columns: ", paste(col.names.repl, collapse=", "),"\n\n",
            "and in the following rows: ", paste(row.names.repl, collapse=", "))
  }

  data.out <- reshape2::dcast(data, proteinname + fullpeptidename_charge  ~ condition + bioreplicate,
                              fun.aggregate=mean, value.var="intensity")
  data.out <- droplevels(data.out)

  return(data.out)
}
