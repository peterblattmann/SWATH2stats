#' Convert data into the format for running a python script.
#'
#' This functions selects the columns suggested to run a python script to change
#' the data from peptide-level to transition-level.
#'
#' The necessary columns are selected and the run column is renamed to filename
#' for the script. The intensities are taken from the column aggr_Peak_Area and
#' therefore the Intensity column is not exported.
#'
#' @param data A data frame containing SWATH data.
#' @param replace.Unimod Option to indicate if Unimod Identifier should be
#'   replaced form ":"" to "_".
#' @return Returns a data frame in the appropriate format to be used by a custom
#'   python script stored in the scripts folder.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data,0.01)
#'  data.pythonscript <- convert4pythonscript(data.filtered.decoy)
#' @export
convert4pythonscript <- function(data, replace.Unimod = TRUE) {

    data <- data[, c("ProteinName", "FullPeptideName", "Charge", "aggr_Fragment_Annotation",
        "aggr_Peak_Area", "RT", "BioReplicate", "Condition", "Run")]

    colnames(data) <- gsub("Run", "filename", colnames(data))

    if (isTRUE(replace.Unimod)) {
        # replace UniMod: to UniMod_
        data$FullPeptideName <- gsub("UniMod:", "UniMod_", data$FullPeptideName)
        data$aggr_Fragment_Annotation <- gsub("UniMod:", "UniMod_", data$aggr_Fragment_Annotation)
    }

    return(data)
}
