#' Take data from the swath2stats format and get it ready for use by MSstats.
#'
#' Though the two tools use very similar formats, some coercion is required to
#' Convert table into the format for MSstats.
#'
#' This functions selects the columns necessary for MSstats and renames them if
#' necessary.
#'
#' The necessary columns are selected and three columns renamed:
#' FullPeptideName -> PeptideSequence
#' Charge -> PrecursorCharge
#' filename -> File
#'
#' @param data  A data frame containing SWATH data.
#' @param replace_values Option to indicate if negative and 0 values should be replaced with NA.
#' @param replace_colnames Option to indicate if column names should be renamed
#'   and columns reduced to the necessary columns for MSstats.
#' @param replace_unimod Option to indicate if Unimod Identifier should be
#'   replaced from ":" to "_".
#' @return Returns a data frame in the appropriate format for MSstats.
#' @references Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek
#'   O. MSstats: an R package for statistical analysis of quantitative mass
#'   spectrometry-based proteomic experiments.Bioinformatics. 2014 Sep
#'   1;30(17):2524-6. doi: 10.1093/bioinformatics/btu305.
#' @author Peter Blattmann
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#'  data.mapDIA <- convert4MSstats(raw)
#' }
#' @export
convert_MSstats <- function(data, replace_values=TRUE,
                                replace_colnames=TRUE, replace_unimod=TRUE) {

  msstats_columns <- c("proteinname", "peptidesequence", "precursorcharge",
                       "fragmention", "productcharge", "isotopelabeltype",
                       "intensity", "bioreplicate", "condition", "run")

  if (isTRUE(replace_colnames)) {
    colnames(data) <- gsub(pattern="fullpeptidename",
                           replacement="peptidesequence",
                           x=colnames(data))
    colnames(data) <- gsub(pattern="^charge$", replacement="precursorcharge",
                           x=colnames(data))
    colnames(data) <- gsub(pattern="align_origfilename", replacement="file",
                           x=colnames(data))
  }
  col_names_missing <- msstats_columns[!(msstats_columns %in% colnames(data))]

  if (length(col_names_missing) > 0) {
      message("One or several columns required by MSstats were not in the data. The columns were created and filled with NAs.\nMissing columns: ",
              paste(unlist(col_names_missing), collapse=", "))

      data[, col_names_missing] <- NA

      if ("isotopelabeltype" %in% col_names_missing) {
        message("isotopelabeltype was filled with light.")
        data[, "isotopelabeltype"] <- "light"
      }

      if ("precursorcharge" %in% col_names_missing) {
          data[, "precursorcharge"] <- gsub(".*_([[:digit:]])$", "\\1", data[, "fragmention"])
      }
    }
    data <- data[, msstats_columns]

    if (isTRUE(replace_values)) {
    # replace negative values to 0 and 0 to NA
      if (sum(data[["intensity"]] < 0, na.rm=TRUE) > 0) {
        data[data[["intensity"]] < 0, "intensity"] <- 0
        warning("Negative intensity values were replaced by NA")
      }

      if (sum(data[["intensity"]] == 0, na.rm=TRUE) > 0) {
        data[data[["intensity"]] %in% 0, "intensity"] <- NA
      warning("Intensity values that were 0, were replaced by NA")
      }
    }

    if(isTRUE(replace_unimod)){
      ## replace UniMod: to UniMod_
      data[["peptidesequence"]] <- gsub("UniMod:", "UniMod_", data[["peptidesequence"]])
      data[["fragmention"]] <- gsub("UniMod:", "UniMod_", data[["fragmention"]])
    }
  return(data)
}
