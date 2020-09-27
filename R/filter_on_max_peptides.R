#' Filter only for the highest intense peptides
#'
#' In order to reduce the data, the data is filtered only for the proteins with
#' the highest intensity peptides.
#'
#' @param data A data frame containing SWATH data with the column names:
#'   ProteinNames, PeptideSequence, PrecursorCharge, Intensity.
#' @param n_peptides Maximum number of highest intense peptides to filter the
#'   data on.
#' @param rm.decoy  Option to remove the decoys during filtering.
#' @return  Returns a data frame of the filtered data.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- filter_mscore_freqobs(data, 0.01,0.8)
#'  data.max <- filter_on_max_peptides(data.filtered, 5)
#' @export
filter_on_max_peptides <- function(data, 
                                   n_peptides, 
                                   rm.decoy = TRUE) {
    data <- unifyProteinGroupLabels(data)
    if (isTRUE(rm.decoy)) {
        data <- removeDecoyProteins(data)
    }

    data <- data.table::data.table(data)
    if (length(grep("ProteinName", colnames(data))) > 0) {
        setnames(data, "ProteinName", "PROTEIN")
    }
    if(length(grep("FullPeptideName", colnames(data))) > 0){
        Peptide_col <- "FullPeptideName"
        setnames(data, "FullPeptideName", "PEPTIDE")
    }
    if(length(grep("PeptideSequence", colnames(data))) > 0){
        Peptide_col <- "PeptideSequence"
        setnames(data, "PeptideSequence", "PEPTIDE")
    }
    
    # data$PEPTIDE <- paste(data$PeptideSequence, data$PrecursorCharge, sep='_')
    data$PEPTIDE <- data$FullPeptideName


    data.peptides <- data[, c("PROTEIN", "PEPTIDE", "Intensity"), with = FALSE]
    data.table::setkeyv(data, cols = c("PROTEIN", "PEPTIDE"))

    Intensity <- NULL
    data.peptides.int <- data.peptides[, sum(Intensity), by = "PROTEIN,PEPTIDE"]
    setnames(data.peptides.int, "V1", "SUM.INTENSITY")

    setkey(data.peptides.int, PROTEIN)
    data.peptides.int <- data.peptides.int[order(data.peptides.int$SUM.INTENSITY,
        decreasing = TRUE), ]

    .SD <- NULL
    peptides.sel <- unique(data.peptides.int[, head(.SD, n_peptides), by = PROTEIN])

    data.filtered <- data.frame(data[PEPTIDE %in% peptides.sel$PEPTIDE, ])


    message("Before filtering: ", "\n", 
            "  Number of proteins: ", length(unique(data$PROTEIN)), "\n", 
            "  Number of peptides: ", length(unique(data$PEPTIDE)), "\n\n", 
            "Percentage of peptides removed: ", 
            round((length(unique(data$PEPTIDE)) - length(unique(data.filtered$PEPTIDE)))/length(unique(data$PEPTIDE)) *
            100, digits = 2), "%", "\n\n", 
            "After filtering: ", "\n", "  Number of proteins: ", legth(unique(data.filtered$PROTEIN)), "\n", 
            "  Number of peptides: ", length(unique(data.filtered$PEPTIDE)),"\n")

    colnames(data.filtered) <- gsub("PROTEIN", "ProteinName", colnames(data.filtered))
    data.filtered <- data.filtered[, -which(colnames(data.filtered) == "PEPTIDE")]

    return(data.filtered)
}
