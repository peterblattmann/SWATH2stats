#' Filter OpenSWATH output table according to mscore and conditions.
#'
#' This function filters the SWATH data according to the m_score value, as well
#' as to the number of occurence in the data (requant) and within a condition
#' (condition).
#'
#' @param data A data frame containing SWATH data.
#' @param mscore Value that defines the mscore threshold according to which the
#'   data will be filtered.
#' @param n.replica Number of measurements within at least one condition that
#'   have to pass the mscore threshold for this transition.
#' @param rm.decoy  Option to drop decoys from the data
#' @param mscore.col Defines the column from which to retrieve the m_score. 
#'    If you use JPP (Rosenberger, Bludau et al. 2017) this can be used to 
#'    select between Protein and transition_group m_score.
#' @return Data which has been filtered.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore(data, 0.01)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01, 0.8)
#'  data.filtered <- filter_mscore_condition(data, 0.01, 3)
#' @export
filter_mscore_condition <- function(data, 
                                    mscore, 
                                    n.replica, 
                                    rm.decoy = TRUE, 
                                    mscore.col = "m_score") {
    .N <- NULL  ## Avoid R CMD check shenanigans

    mscore.col <- JPP_update(data, mscore.col)

    # a column with unique identifiers for each Precursor (e.g. ADFSDF 3) 
    # is generated
    data$Peptide_Charge <- paste(data$FullPeptideName, data$Charge)
    ## a column with unique identifiers for each Precursor and Condition is
    ## genrated ADFSDF 3 Condition 1
    data$Peptide_Charge_Condition <- paste(data$Peptide_Charge, data$Condition)

    ## decoys are removed if present
    if (sum(colnames(data) == "decoy") == 1 & isTRUE(rm.decoy)) {
        data <- data[data$decoy == 0, ]
    }

    ## only data that is below the indicated mscore is selected and then only 
    ## the unique data selected data.filtered <- subset(data, m_score <= mscore)
    data.filtered <- data[data[, mscore.col] <= mscore, ]
    data.filtered <- unique(data.filtered[, c("Peptide_Charge", 
                                              "Peptide_Charge_Condition",
                                              "aggr_Peak_Area")])
    data.filtered <- data.table::data.table(data.filtered)

    data.table::setkeyv(data.filtered,
                        cols = c("Peptide_Charge", "Peptide_Charge_Condition", 
                                 "aggr_Peak_Area"))
    ## number of occurences of Precursor per Condition is calculated
    data.n <- data.filtered[, .N, by = "Peptide_Charge,Peptide_Charge_Condition"]

    ## only precursors that are present in more that n.replica are selected
    precursor.filtered <- data.n[data.n$N >= n.replica]
    precursor.filtered <- data.frame(
      "Peptide_Charge" = unique(precursor.filtered$Peptide_Charge))

    message("Fraction of peptides selected: ",
            signif(length(unique(precursor.filtered$Peptide_Charge)) /
                   length(unique(data$Peptide_Charge)),
                   digits = 2))

    # only data that is present in the precursor.filtered list is selected
    data.filtered <- merge(data, precursor.filtered,
                           by.x = "Peptide_Charge", by.y = "Peptide_Charge")

    message("Dimension difference: ", paste(dim(data) - dim(data.filtered), 
                                            collapse = ", "))
    return(data.filtered)
}
