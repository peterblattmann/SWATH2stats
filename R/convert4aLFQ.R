#' Convert table into the format expected by aLFQ.
#'
#' This function selects the columns necessary for the aLFQ R package.
#'
#' @param data A data frame containing the SWATH data in transition-level
#'   format
#' @param annotation Option to indicate if the data has been annotated, i.e. if
#'   the columns Condition, Replicate, Run are present. If option
#'   is set to true it will write a new run_id as a string of the
#'   combination of these three columns.
#' @param check_transitions Option if number of transitions should be checked.
#'   As input only transition-level data should be used and
#'   therefore this is checked. However, this makes the function
#'   slow and herewith be omitted.
#' @return Returns a data frame in the appropriate format for aLFQ.
#' @author Peter Blattmann
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- SWATH2stats::sample_annotation(OpenSWATH_data, Study_design, verbose=TRUE)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#'  data.aLFQ <- convert4aLFQ(raw)
#'  }
#' @importFrom stats aggregate median
#' @export

convert4aLFQ <- function(data, 
                         annotation = TRUE, 
                         check_transitions = TRUE) {

    if (annotation == TRUE) {
        data <- data[, c("ProteinName", "PeptideSequence", "FragmentIon", 
                         "NakedSequence", "PrecursorCharge", "Intensity", 
                         "Condition", "BioReplicate", "Run")]
        data$run_id <- paste(data$Condition, data$BioReplicate, data$Run, sep = "_")
    }

    colnames(data) <- gsub("ProteinName", "protein_id", colnames(data))
    colnames(data) <- gsub("PrecursorCharge", "precursor_charge", colnames(data))
    colnames(data) <- gsub("FragmentIon", "transition_id", colnames(data))
    colnames(data) <- gsub("Intensity", "transition_intensity", colnames(data))
    colnames(data) <- gsub("PeptideSequence", "peptide_id", colnames(data))
    colnames(data) <- gsub("NakedSequence", "peptide_sequence", colnames(data))

    data$transition_id <- gsub("_run[[:digit:]]*$", "", data$transition_id)
    data$transition_id <- paste(data$peptide_sequence, data$transition_id, sep = " ")

    data$concentration <- "?"

    data <- data[, c("run_id", "protein_id", "peptide_id", "transition_id", 
                     "peptide_sequence", "precursor_charge", 
                     "transition_intensity", "concentration")]
    # check transitions
    data$protein_id <- factor(data$protein_id)
    data$peptide_id <- factor(data$peptide_id)
    data$transition_id <- factor(data$transition_id)

    if (check_transitions) {
        message("Checking the integrity of the transitions takes a lot of time. 
                To speed up consider changing the option.")
        data.agg <- aggregate(data[, c("transition_id")], by = list(data$peptide_id,
            data$run_id), length)
        if (median(data.agg$x) == 1) {
            warning("The aLFQ package should only be used with transition-level data.
              The data only contains one transition per peptide.")
        }
    }

    # convert back to character vector
    data$protein_id <- as.character(data$protein_id)
    data$peptide_id <- as.character(data$peptide_id)
    data$transition_id <- as.character(data$transition_id)

    return(data)
}
