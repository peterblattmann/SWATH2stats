#' Reduce columns of OpenSWATH data
#'
#' This function selects the columns from the standard OpenSWATH output to
#' column needed for MSstats, aLFQ and mapDIA.
#'
#' @note A basic set of columns are defined in the function and are used if no column
#'   names are indicated.
#' @note The column.names can be omitted and then the following columns are
#'   selected that are needed for MSstats and mapDIA analysis: ProteinName,
#'   FullPeptideName, Sequence, Charge, aggr_Fragment_Annotation,
#'   aggr_Peak_Area, filename, m_score, decoy, Intensity, RT.
#'   This function should be ommitted if the data is analyzed afterwards with
#'   the aLFQ or imsbInfer package that needs further columns.
#'
#' @param data  A data frame containing SWATH data.
#' @param column.names A vector of column names that can be selected.
#' @return Returns a data frame with the selected columns.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered <- reduce_OpenSWATH_output(data)
#' @export
reduce_OpenSWATH_output <- function(data, 
                                    column.names = NULL) {
    if (is.null(column.names)) {
        if ("transition_group_id_m_score" %in% colnames(data)) {
            column.names <- c("ProteinName", "FullPeptideName", "Sequence", 
                              "Charge","aggr_Fragment_Annotation", 
                              "aggr_Peak_Area", "filename", 
                              "transition_group_id_m_score", 
                              "ProteinName_m_score", "decoy", "Intensity", 
                              "RT", "run_id", "transition_group_id")
        }
        if ("m_score" %in% colnames(data)) {
            column.names <- c("ProteinName", "FullPeptideName", "Sequence", 
                              "Charge", "aggr_Fragment_Annotation", 
                              "aggr_Peak_Area", "filename", "m_score",
                              "decoy", "Intensity", "RT", "run_id", 
                              "transition_group_id")
        }

    }
    if (length(column.names) > length(column.names[column.names %in% colnames(data)])) {
        col.names.missing <- column.names[!column.names %in% colnames(data)]
        warning("These columns are missing from the data:", paste(unlist(col.names.missing),
            collapse = ", "))

    }
    # Keep only required columns for MSStats and mapDIA
    if (length(column.names) == length(column.names[column.names %in% colnames(data)])) {
        data.filtered <- data[, column.names]
        return(data.filtered)
    }
}
