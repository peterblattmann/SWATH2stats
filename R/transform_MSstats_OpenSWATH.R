#'  Transforms column names to OpenSWATH column names
#'
#' This functions transforms the column names from a data frame in MSstats
#' format to a data frame with column names used by the OpenSWATH output. The
#' original table needs to contain at least the 10 columns defined by MSstats:
#' ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
#' IsotopeLabelType, Condition, BioReplicate, Run, Intensity.)
#'
#' @param data A data frame containing the SWATH data in the MSstats format
#' @return The data frame in the appropriate format.
#' @references Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek
#'   O. MSstats: an R package for statistical analysis of quantitative mass
#'   spectrometry-based proteomic experiments.Bioinformatics. 2014 Sep
#'   1;30(17):2524-6. doi: 10.1093/bioinformatics/btu305.
#' @author Peter Blattmann
#' @examples
#' \dontrun{
#'  data("MSstats_data", package="SWATH2stats")
#'  transform_MSstats_OpenSWATH(MSstats_data)
#' }
#' @export
transform_MSstats_OpenSWATH <- function(data) {
  colnames(data) <- gsub("\\.", "", colnames(data))
  colnames(data) <- tolower(colnames(data))

  missing_columns <- setdiff(c("proteinname", "peptidesequence", "precursorcharge",
                               "fragmention", "productcharge", "isotopelabeltype",
                               "condition", "bioreplicate", "run", "intensity"), colnames(data))

  if (length(missing_columns > 0)) {
    stop("The data frame doesn't contain all required columns. The following columns are missing:",
         paste(missing_columns, collapse=", "))
  }

  query_columns <- c("proteinname", "peptidesequence", "precursorcharge", "fragmention",
                     "productcharge", "isotopelabeltype", "condition", "bioreplicate",
                     "run", "intensity")
  add.colnames <- colnames(data)[!(colnames(data) %in% query_column)]

  data <- data[, c(query_columns, add.colnames)]
  new_columns <- gsub(pattern="peptidesequence", replacement="fullpeptidename", x=query_columns)
  colnames(data)[seq_len(10)] <- new_columns

  if (length(add.colnames) > 0) {
    message("Additional columns present in the data: ", paste(add.colnames, collapse=", "))
  }

  if (!("mscore" %in% colnames(data))) {
    message("No column 'mscore' present in the data. This column is required for functions estimating FDR or filtering the data.")
  }
  return(data)
}
