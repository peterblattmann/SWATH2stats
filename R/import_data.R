#' Transforms the column names from a data frame to the required format.
#'
#' This functions transforms the column names from a data frame from another
#' format to a data frame with column names used by the OpenSWATH output and
#' required for these functions. During executing of the function the
#' corresponding columns for each column in the data need to be selected. For
#' columns that do not corresond to a certain column 'not applicable' needs to
#' be selected and the column names are not changed.
#'
#' @note List of column names of the OpenSWATH data:
#' ProteinName: Unique identifier for protein or proteingroup that the peptide
#'   maps to. Proteotypic peptides should be indicated by 1/ in order to be
#'   recognized as such by the function filter_proteotypic_peptides.
#' FullPeptideName: Unique identifier for the peptide.
#' Charge: Charge of the peptide precursor ion quantified.
#' Sequence: Naked peptide sequence without modifications.
#' aggr_Fragment_Annotation: aggregated annotation for the different Fragments
#'   quantified for this peptide. In the OpenSWATH results the different
#'   annotation in OpenSWATH are concatenated by a semicolon.
#' aggr_Peak_Area: aggregated Intensity values for the different Fragments
#'   quantified for this peptide. In the OpenSWATH results the aggregated Peak
#'   Area intensities are concatenated by a semicolon.
#' transition_group_id: A unique identifier for each transition group used.
#' decoy: Indicating with 1 or 0 if this transition group is a decoy.
#' m_score: Column containing the score that is used to estimate FDR or
#'   filter. M-score values of identified peak groups are equivalent to a
#'   q-value and thus typically are smaller than 0.01, depending on the
#'   confidence of identification (the lower the m-score, the higher the
#'   confidence).
#' Column containing the score that is used to estimate FDR or filter.
#' RT: Column containing the retention time of the quantified peak.
#' filename: Column containing the filename or a unique identifier for each
#'   injection.
#' Intensity: column containing the intensity value for each quantified peptide.
#' Columns needed for FDR estimation and filtering functions: ProteinName,
#'   FullPeptideName, transition_group_id, decoy, m_score
#' Columns needed for conversion to transition-level format (needed for MSStats
#'   and mapDIA input): aggr_Fragment_Annotation, aggr_Peak_Are
#'
#' @param data  A data frame containing the SWATH-MS data (one line per peptide
#'   precursor quantified) but with different column names.
#' @return Returns the data frame in the appropriate format.
#' @author Peter Blattmann
#' @examples
#'  data('Spyogenes', package = 'SWATH2stats')
#'  head(data)
#'  str(data)
#' @export
import_data <- function(data) {
  colnames_OpenSWATH <- c("proteinname", "fullpeptidename", "charge", "sequence",
                          "aggr_fragment_annotation", "aggr_peak_area",
                          "transition_group_id", "decoy", "m_score",
                          "rt", "align_origfilename", "intensity", "not applicable")

  colnames(data) <- tolower(colnames(data))
  message("When reading in data, please specify to which column in the OpenSWATH data they correspond.
For columns that do not correspond to any OpenSWATH column, choose \"not applicable\".
Explanation of the OpenSWATH columns can be found in the manual page.")

  # Dialogue to map columns
  for (i in colnames(data)) {
    value <- select.list(
      colnames_OpenSWATH,
      title=paste0("Select a column from the OpenSWATH output that the column ", i, "corresponds to."))
    if (value != "not applicable") {
      colnames(data) <- gsub(pattern=i, replacement=value, x=colnames(data))
      message("Column name ", i, " was replaced with ", value, ".")
    }
  }

  # add NA to colnames that were not mapped
  add.colnames <- colnames(data)[!(colnames(data) %in% colnames_OpenSWATH)]
  if (length(add.colnames > 0)) {
    for(i in add.colnames) {
      data[, i] <- NA
    }
    warning("Not all columns required within SWATH2stats were mapped: The columns ",
            paste(add.colnames, collapse =", "), " were added with NA as value.")
  }

  return(data)
}
