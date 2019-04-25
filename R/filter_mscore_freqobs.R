utils::globalVariables(c("Peptide_Charge", ".N"))

#' Filter openSWATH output table according to mscore.
#'
#' This function filters the SWATH data according to the m_score value, as well
#' as to the number of occurence in the data (requant) and within a condition
#' (condition)
#'
#' @param data  A data frame containing SWATH data.
#' @param mscore  Value that defines the mscore threshold according to which the
#'   data will be filtered.
#' @param percentage  Percentage in which replicas the transition has to reach
#'   the mscore threshold. ??
#' @param rm.decoy  Option to remove the decoys during filtering.
#' @param file_column  Column in the (raw)data in which the filenames should be found.
#' @return Returns a data frame with the filtered data.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore(data, 0.01)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01, 0.8)
#'  data.filtered <- filter_mscore_condition(data, 0.01, 3)
#' @export

filter_mscore_freqobs <- function(data, mscore, percentage = NULL, rm.decoy = TRUE,
    mscore.col = "m_score") {

    mscore.col <- JPP_update(data, mscore.col)

    data$Peptide_Charge <- paste(data$FullPeptideName, data$Charge)

    if (sum(colnames(data) == "decoy") == 1 & isTRUE(rm.decoy)) {
        # data <- subset(data, decoy == 0)
        data <- data[data$decoy == 0, ]
    }

    # data.filtered <- subset(data, m_score <= mscore)
    data.filtered <- data[data[, mscore.col] <= mscore, ]
    data.filtered <- data.table(data.filtered)

    data.filtered <- data.filtered[, c("Peptide_Charge", "aggr_Peak_Area"), with = FALSE]
    setkey(data.filtered, Peptide_Charge)
    data.n <- data.filtered[, .N, by = "Peptide_Charge"]

    if (is.null(percentage)) {
        percentage <- 0
    }

    threshold <- nlevels(factor(data$filename)) * percentage
    message("Treshold, peptides need to have been quantified in more conditions than: ",
        threshold)

    peptides.filtered <- data.n[data.n$N >= threshold]
    peptides.filtered <- data.frame(Peptides_Charge = peptides.filtered$Peptide_Charge)

    message("Fraction of peptides selected: ", signif(length(unique(peptides.filtered$Peptides_Charge))/length(unique(data$Peptide_Charge)),
        digits = 2))

    data.filtered <- merge(data, peptides.filtered, by.x = "Peptide_Charge", by.y = "Peptides_Charge")

    message("Dimension difference: ", paste(dim(data) - dim(data.filtered), collapse = ", "))

    return(data.filtered)
}
