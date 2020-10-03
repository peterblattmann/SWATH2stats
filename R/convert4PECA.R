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
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  data.PECA <- convert4PECA(data.filtered.decoy)
#' }
#' @importFrom reshape2 dcast
#' @export

convert4PECA <- function(data) {

    columns <- c("ProteinName", "FullPeptideName", "Charge", "Intensity", 
                 "BioReplicate", "Condition", "RT")

    missing_columns <- columns[!(columns %in% colnames(data))]
    if (length(missing_columns) > 0) {
        warning("One or several columns required by mapDIA were not in the data and filled with NAs.
            Missing columns: ",
            paste(unlist(missing_columns), collapse = ", "))

        data[, missing_columns] <- NA
    }

    data$FullPeptideName_Charge <- paste(data$FullPeptideName, data$Charge, sep = "_")

    # # replace UniMod: to UniMod_ data$PeptideSequence <- gsub(':', '_',
    # data$PeptideSequence) data$FragmentIon <- gsub(':', '_', data$FragmentIon)

    data.wide <- dcast(data, ProteinName + FullPeptideName_Charge ~ Condition +
        BioReplicate, fun.aggregate = length, value.var = "Intensity")

    if (sum(data.wide[, 4:dim(data.wide)[2]] > 1) > 0) {
        duplicates.cols <- apply(data.wide[, 4:dim(data.wide)[2]],
            2, function(x) sum(x > 1))
        duplicates.cols <- names(duplicates.cols[duplicates.cols > 0])
        duplicates.rows <- apply(data.wide[, 4:dim(data.wide)[2]],
            1, function(x) sum(x > 1))
        duplicates.rows <- data.wide[duplicates.rows > 0, "FullPeptideName_Charge"]

        warning("Data contains several intensity values per condition\n\n", "in the following columns: ",
            paste(duplicates.cols, collapse = ", "), "\n\n", "and in the following rows: ",
            paste(duplicates.rows, collapse = ", "))
    }

    data.out <- dcast(data, ProteinName + FullPeptideName_Charge ~ Condition + BioReplicate,
        fun.aggregate = mean, value.var = "Intensity")

    data.out <- droplevels(data.out)

    return(data.out)
}
