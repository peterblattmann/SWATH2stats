#' Transforms the SWATH data from a peptide- to a transition-level table.
#'
#' If the SWATH data should be analyzed on transition-level the data needs to be
#' tranformed from peptide-level table to a transition-level table (one row per
#' transition instead of one row per peptide). The columns
#' "aggr_Fragment_Annotation" and "aggr_Peak_Area" are disaggregated into the
#' new columns "FragmentIon" and "Intensity".
#' The following columns are renamed if they exist: FullPeptideName ->
#' PeptideSequence, Charge -> PrecursorCharge, Area -> Intensity, Fragment ->
#' FragmentIon, Sequence -> NakedSequence.
#'
#' @param data A data frame containing SWATH data.
#' @param all.columns Option that all columns are processed. Otherwise only the
#'   typical columns needed for downstream analysis are processed.
#' @return Returns a data frame containing the SWATH data in a transition-level
#'   table.
#' @author Peter Blattmann
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#' }
#' @export
disaggregate <- function(data, all.columns=FALSE) {

  # sanity test on the number of transitions per precursor
  n.transitions <- lapply(as.character(data[["aggr_fragment_annotation"]]),
                          function(x)strsplit(x, ";"))
  n.transitions2 <- unlist(lapply(n.transitions, function(x) length(unlist(x))))

  n.transitions3 <- lapply(as.character(data[["aggr_peak_area"]]), function(x)strsplit(x, ";"))
  n.transitions4 <- unlist(lapply(n.transitions3, function(x) length(unlist(x))))

  if (sum(n.transitions2 != n.transitions4) > 0) {
    stop(paste("The number of transitions annotated and measured do not match in the following transitions:\n",
               paste(unlist(n.transitions[n.transitions2 != n.transitions4]), collapse=", ")))
  }

  # test if always the same number of transitions per precursor were used
  if (min(n.transitions2) == max(n.transitions2)) {
    message("The library contains", max(n.transitions2), "transitions per precursor.
The data table was transformed into a table containing one row per transition.")
  }

  if (min(n.transitions2) != max(n.transitions2)) {
    message("The library contains between ", min(n.transitions2)," and ", max(n.transitions2),
            " transitions per precursor.
The data table was transformed into a table containing one row per transition.")
  }

  data.new <- cbind(data, reshape2::colsplit(data[["aggr_fragment_annotation"]], ";",
                                             paste("split_fragannot_",
                                                   seq_len(max(n.transitions2)), sep="")),
                    reshape2::colsplit(data[["aggr_peak_area"]], ";",
                                       paste("split_peakarea_",
                                             seq_len(max(n.transitions2)), sep="")), stringsAsFactors=FALSE)

  data.new.m <- reshape2::melt(data.new, id.vars=grep("split_fragannot",
                                                      colnames(data.new), invert=TRUE),
                               measure.vars=grep("split_fragannot", colnames(data.new)),
                               variable.name="fragannot_n", value.name="fragment")
  data.new.m2 <- reshape2::melt(data.new.m, id.vars=grep("split_peakarea",
                                                         colnames(data.new.m), invert=TRUE),
                                measure.vars=grep("split_peakarea", colnames(data.new.m)),
                                variable.name="area_n", value.name="area")

  # added because it didn't name the variables in a later trial
  if (sum(colnames(data.new.m2) %in% c("fragannot_n")) == 0) {
    l <- length(colnames(data.new.m2))
    colnames(data.new.m2)[l - 3] <- "fragannot_n"
    colnames(data.new.m2)[l - 2] <- "fragment"
    colnames(data.new.m2)[l - 1] <- "area_n"
    colnames(data.new.m2)[l] <- "area"
  }

  data.new.m3 <- data.new.m2[gsub("split_fragannot_", "",
                                  data.new.m2[, "fragannot_n"]) == gsub("split_peakarea_", "",
                                                                       data.new.m2[,"area_n"]), ]
  cols <- colnames(data.new.m3)
  if (isTRUE(all.columns)) {
    data.new.merged <- data.new.m3
  } else {
    cols <- colnames(data.new.m3)[colnames(data.new.m3) %in%
                                  c("proteinname", "fullpeptidename", "peptidesequence", "sequence",
                                    "charge", "precursorcharge","fragment", "fragmention",
                                    "area", "condition", "bioreplicate", "run", "rt")]
    data.new.merged <- data.new.m3[, cols]
  }

  colnames(data.new.merged) <- gsub("fullpeptidename", "peptidesequence", colnames(data.new.merged))
  colnames(data.new.merged) <- gsub("^charge$", "precursorcharge", colnames(data.new.merged))

  colnames(data.new.merged) <- gsub("area", "intensity", colnames(data.new.merged))
  colnames(data.new.merged) <- gsub("fragment", "fragmention", colnames(data.new.merged))

  if ("sequence" %in% cols) {
    colnames(data.new.merged) <- gsub("^sequence$", "nakedsequence", colnames(data.new.merged))
  }

  if (sum(is.na(data.new.merged[["intensity"]])) > 0) {
    .ids <- !is.na(data.new.merged[["intensity"]])
    message(length(data.new.merged[["intensity"]]) - sum(.ids),
            " row(s) was(were) removed because they did not contain data due to different number of transitions per precursor.")
    data.new.merged <- data.new.merged[.ids, ]
  }

  return(data.new.merged)
}
