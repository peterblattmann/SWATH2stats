#' Filter openSWATH output for proteins that are identified by a minimum of n independent peptides.
#'
#' This function removes entries mapping to proteins that are identified by less
#' than n_peptides.
#' Removing single-hit proteins from an analysis can significantly increase the
#' sensitivity under strict protein fdr criteria, as evaluated by
#' e.g. assess_fdr_overall.
#'
#' @param data Data table that is produced by the OpenSWATH/iPortal workflow.
#' @param n_peptides Number of minimal number of peptide IDs associated with a
#'   protein ID in order to be kept in the dataset.
#' @param protein_col Column with protein identifiers. Default: ProteinName
#' @param peptide_col Column with peptide identifiers. Default: Peptide.Sequence or FullPeptideName
#' @param rm.decoy  Option to remove the decoys during filtering.
#' @return Returns the filtered data frame with only peptides that map to
#'   proteins with >= n_peptides peptides.
#' @author Moritz Heusel, Peter Blattmann
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01,0.8)
#'  data.max <- filter_on_max_peptides(data.filtered, 5)
#'  data.min.max <- filter_on_min_peptides(data.max, 3)
#'  }
#' @export
filter_on_min_peptides <- function(data, 
                                   n_peptides, 
                                   protein_col = "ProteinName",
                                   peptide_col = c("Peptide.Sequence", "FullPeptideName"),
                                   rm.decoy = TRUE) {
    data <- unifyProteinGroupLabels(data)
    if (isTRUE(rm.decoy)) {
        data <- removeDecoyProteins(data)
    }
    
    # select valid columns
    columns <- validate_columns(data, list(Protein = protein_col,
                                        Peptide = peptide_col))
    
    # rename columns
    colnames(data) <- gsub(columns[["Protein"]], "PROTEIN", colnames(data))
    colnames(data) <- gsub(columns[["Peptide"]], "PEPTIDE", colnames(data))
    
    data.prot.pep <- unique(data[, c("PROTEIN", "PEPTIDE")])

    data.prot.pep.n <- tapply(data.prot.pep$PEPTIDE, 
                              data.prot.pep$PROTEIN,
                              length)

    prot.pep.names <- names(data.prot.pep.n[data.prot.pep.n >= n_peptides & 
                                                !is.na(data.prot.pep.n)])
    data.filtered <- data[data$PROTEIN %in% prot.pep.names, ]

    message("Before filtering: ", "\n", 
            "  Number of proteins: ", length(unique(data$PROTEIN)), "\n", 
            "  Number of peptides: ", length(unique(data$PEPTIDE)), "\n\n",
            "Percentage of peptides removed: ", round((length(unique(data$PEPTIDE)) -
            length(unique(data.filtered$PEPTIDE)))/length(unique(data$PEPTIDE)) *
            100, digits = 2), "%", "\n\n", 
            "After filtering: ", "\n", "  Number of proteins: ",
            length(unique(data.filtered$PROTEIN)), "\n", 
            "  Number of peptides: ", length(unique(data.filtered$PEPTIDE)), "\n")
   
    # rename columns to original names
    colnames(data.filtered) <- gsub("PROTEIN", columns[["Protein"]], colnames(data.filtered))
    colnames(data.filtered) <- gsub("PEPTIDE", columns[["Peptide"]], colnames(data.filtered))
    
    return(data.filtered)
}
