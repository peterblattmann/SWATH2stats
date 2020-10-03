#' Filter OpenSWATH output table according to mscore and conditions.
#'
#' This function filters the SWATH data according to the m_score value, as well
#' as to the number of occurence in the data (requant) and within a condition
#' (condition).
#'
#' @param data A data frame containing SWATH data.
#' @param mscore Value that defines the mscore threshold according to which the
#'   data will be filtered.
#' @param n_replica Number of measurements within at least one condition that
#'   have to pass the mscore threshold for this transition.
#' @param peptide_col Column with peptide identifiers. Default: Peptide.Sequence or FullPeptideName
#' @param charge_col Column with peptide charge. Default: Charge
#' @param condition_col Column with conditions. Default: Condition
#' @param rm.decoy  Option to drop decoys from the data
#' @param mscore.col Defines the column from which to retrieve the m_score. 
#'    If you use JPP (Rosenberger, Bludau et al. 2017) this can be used to 
#'    select between Protein and transition_group m_score.
#' @return Data which has been filtered.
#' @author Peter Blattmann
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore(data, 0.01)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01, 0.8)
#'  data.filtered <- filter_mscore_condition(data, 0.01, 3)
#'  }
#' @importFrom data.table setkeyv data.table
#' @export
filter_mscore_condition <- function(data, 
                                    mscore, 
                                    n_replica, 
                                    peptide_col = c("Peptide.Sequence", "FullPeptideName"),
                                    charge_col = "Charge",
                                    condition_col = "Condition",
                                    rm.decoy = TRUE, 
                                    mscore.col = "m_score") {
    .N <- NULL  ## Avoid R CMD check shenanigans

    mscore.col <- JPP_update(data, mscore.col)
   
     ## decoys are removed if present
    if (sum(colnames(data) == "decoy") == 1 & isTRUE(rm.decoy)) {
      data <- data[data$decoy == 0, ]
    }
    
    #Select valid columns
    columns <- validate_columns(data, list(Peptide = peptide_col,
                                        Charge = charge_col,
                                        Condition = condition_col,
                                        Score = mscore.col))
    
    
    # a unique identifiers for each Precursor (e.g. ADF_3) is generated
    data$Peptide_Charge <- paste(data[,columns[["Peptide"]]], data[,columns[["Charge"]]])
    ## a unique identifiers for each Precursor and Condition is generated 
    data$Peptide_Charge_Condition <- paste(data$Peptide_Charge, data[,columns[["Condition"]]])

    ## only data that is below the indicated mscore is selected and then only 
    ## the unique data selected data.filtered <- subset(data, m_score <= mscore)
    data.filtered <- data[data[, columns[["Score"]]] <= mscore, ]
    data.filtered <- unique(data.filtered[, c("Peptide_Charge", 
                                              "Peptide_Charge_Condition",
                                              "aggr_Peak_Area")])
    data.filtered <- data.table::data.table(data.filtered)

    # aggr_Peak area is only used to have a unique column
    data.table::setkeyv(data.filtered,
                        cols = c("Peptide_Charge", "Peptide_Charge_Condition", 
                                 "aggr_Peak_Area"))
    
    ## number of occurrences of Precursor per Condition is calculated
    data.n <- data.filtered[, .N, by = "Peptide_Charge,Peptide_Charge_Condition"]

    ## only precursors that are present in more that n_replica are selected
    precursor.filtered <- data.n[data.n$N >= n_replica]
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
