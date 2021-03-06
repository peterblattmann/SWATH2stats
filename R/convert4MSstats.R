#' Get data ready for use by MSstats.
#'
#' Though SWATH2stats uses very similar format as MSstats, some coercion is 
#' required to convert the data into the format for MSstats.
#'
#' This functions selects the columns necessary for MSstats and renames them if
#' necessary.
#'
#' The necessary columns are selected and three columns renamed:
#' FullPeptideName -> PeptideSequence
#' Charge -> PrecursorCharge
#' filename -> File
#'
#' @param data  A data frame containing SWATH data.
#' @param replace_values Option to indicate if negative and 0 values should be replaced with NA.
#' @param replace_colnames Option to indicate if column names should be renamed
#'   and columns reduced to the necessary columns for MSstats.
#' @param replace_unimod Option to indicate if Unimod Identifier should be
#'   replaced from ":" to "_".
#' @return Returns a data frame in the appropriate format for MSstats.
#' @references Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek
#'   O. MSstats: an R package for statistical analysis of quantitative mass
#'   spectrometry-based proteomic experiments.Bioinformatics. 2014 Sep
#'   1;30(17):2524-6. doi: 10.1093/bioinformatics/btu305.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#'  data.mapDIA <- convert4MSstats(raw)
#' @export

convert4MSstats <- function(data, 
                            replace_values = TRUE, 
                            replace_colnames = TRUE,
                            replace_unimod = TRUE) {

    MSstats_columns <- c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                         "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                         "Intensity", "BioReplicate", "Condition", "Run")

    if (isTRUE(replace_colnames)) {
        colnames(data) <- gsub("FullPeptideName", "PeptideSequence", colnames(data))
        colnames(data) <- gsub("^Charge$", "PrecursorCharge", colnames(data))
        colnames(data) <- gsub("filename", "File", colnames(data))
    }
    missing_columns <- MSstats_columns[!(MSstats_columns %in% colnames(data))]
    if (length(missing_columns) > 0) {
        message("One or several columns required by MSstats were not in the data. 
                The columns were created and filled with NAs.\nMissing columns: ",
            paste(unlist(missing_columns), collapse = ", "))

        if ("IsotopeLabelType" %in% missing_columns) {
            message("IsotopeLabelType was filled with light.")
        }

        data[, missing_columns] <- NA

        if ("IsotopeLabelType" %in% missing_columns) {
            data[, "IsotopeLabelType"] <- "light"
        }
        if ("PrecursorCharge" %in% missing_columns) {
            data[, "PrecursorCharge"] <- gsub(".*_([[:digit:]])$", "\\1", 
                                              data[,"FragmentIon"])
        }
    }
    data <- data[, MSstats_columns]

    if (isTRUE(replace_values)) {
        # replace negative values to 0 and 0 to NA
        if (sum(data$Intensity < 0, na.rm = TRUE) > 0) {
            data[data$Intensity < 0, "Intensity"] <- 0
            warning("Negative intensity values were replaced by NA")
        }
        if (sum(data$Intensity == 0, na.rm = TRUE) > 0) {
            data[data$Intensity %in% 0, "Intensity"] <- NA
            warning("Intensity values that were 0, were replaced by NA")
        }
    }

    if (isTRUE(replace_unimod)) {
        # replace UniMod: to UniMod_
        data$PeptideSequence <- gsub("UniMod:", "UniMod_", data$PeptideSequence)
        data$FragmentIon <- gsub("UniMod:", "UniMod_", data$FragmentIon)
    }
    return(data)
}
