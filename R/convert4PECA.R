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
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.PECA <- convert_PECA(data.filtered.decoy)
#' @export

convert4PECA <- function(data) {

    columns <- c("ProteinName", "FullPeptideName", "Charge", "Intensity", 
                 "BioReplicate", "Condition", "RT")

    col.names.missing <- columns[!(columns %in% colnames(data))]
    if (length(col.names.missing) > 0) {
        warning("One or several columns required by mapDIA were not in the data and filled with NAs.
            Missing columns: ",
            paste(unlist(col.names.missing), collapse = ", "))

        data[, col.names.missing] <- NA
    }

    data$FullPeptideName_Charge <- paste(data$FullPeptideName, data$Charge, sep = "_")

    # # replace UniMod: to UniMod_ data$PeptideSequence <- gsub(':', '_',
    # data$PeptideSequence) data$FragmentIon <- gsub(':', '_', data$FragmentIon)

    data.red.wide.test <- dcast(data, ProteinName + FullPeptideName_Charge ~ Condition +
        BioReplicate, fun.aggregate = length, value.var = "Intensity")

    if (sum(data.red.wide.test[, 4:dim(data.red.wide.test)[2]] > 1) > 0) {
        col.names.repl <- apply(data.red.wide.test[, 4:dim(data.red.wide.test)[2]],
            2, function(x) sum(x > 1))
        col.names.repl <- names(col.names.repl[col.names.repl > 0])
        row.names.repl <- apply(data.red.wide.test[, 4:dim(data.red.wide.test)[2]],
            1, function(x) sum(x > 1))
        row.names.repl <- data.red.wide.test[row.names.repl > 0, "FullPeptideName_Charge"]

        warning("Data contains several intensity values per condition\n\n", "in the following columns: ",
            paste(col.names.repl, collapse = ", "), "\n\n", "and in the following rows: ",
            paste(row.names.repl, collapse = ", "))
    }

    data.out <- dcast(data, ProteinName + FullPeptideName_Charge ~ Condition + BioReplicate,
        fun.aggregate = mean, value.var = "Intensity")

    data.out <- droplevels(data.out)

    return(data.out)
}
