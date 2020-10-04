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
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#'  data.mapDIA <- convert4mapDIA(raw, RT=TRUE)
#'  }
#' @importFrom reshape2 dcast
#' @importFrom stats aggregate
#' @export

convert4mapDIA <- function(data, RT = FALSE) {

    mapDIA_columns <- c("ProteinName", "PeptideSequence", "FragmentIon", 
                      "Intensity", "BioReplicate", "Condition", "RT")

    missing_columns <- mapDIA_columns[!(mapDIA_columns %in% colnames(data))]
    if (length(missing_columns) > 0) {
        warning("One or several columns required by mapDIA were not in the data and filled with NAs.
        Missing columns: ",  paste(unlist(missing_columns), collapse = ", "))
        data[, missing_columns] <- NA
    }

    # replace UniMod: to UniMod_
    data$PeptideSequence <- gsub(":", "_", data$PeptideSequence)
    data$FragmentIon <- gsub(":", "_", data$FragmentIon)

    data.wide.test <- dcast(data, ProteinName + PeptideSequence + FragmentIon ~
        Condition + BioReplicate, fun.aggregate = length, value.var = "Intensity")

    if (sum(data.wide.test[, 4:dim(data.wide.test)[2]] > 1) > 0) {
        col.names.repl <- apply(data.wide.test[, 4:dim(data.wide.test)[2]],
            2, function(x) sum(x > 1))
        col.names.repl <- names(col.names.repl[col.names.repl > 0])
        row.names.repl <- apply(data.wide.test[, 4:dim(data.wide.test)[2]],
            1, function(x) sum(x > 1))
        row.names.repl <- data.wide.test[row.names.repl > 0, "FragmentIon"]

        warning("Data contains several intensity values per condition\n\n", "in the following columns: ",
            paste(col.names.repl, collapse = ", "), "\n\n", "and in the following rows: ",
            paste(row.names.repl, collapse = ", "))
    }

    data.wide <- dcast(data, ProteinName + PeptideSequence + FragmentIon ~ Condition +
        BioReplicate, fun.aggregate = mean, value.var = "Intensity")
    if (RT) {
        RTs <- unique(aggregate(data[, c(which(colnames(data) == "RT"))], by = list(data$FragmentIon),
            FUN = mean, na.rm = TRUE))
        colnames(RTs) <- c("FragmentIon", "RT")
        RTs$RT <- RTs$RT/60
        data.wide <- merge(data.wide, RTs, by = "FragmentIon", all.x = TRUE)
        data.wide <- data.wide[, c("ProteinName", "PeptideSequence", "FragmentIon",
            colnames(data.wide)[!(colnames(data.wide) %in% c("ProteinName",
                "PeptideSequence", "FragmentIon", "RT"))], "RT")]
    }
    return(data.wide)
}
