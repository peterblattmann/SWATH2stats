#' Transforms the SWATH data from a peptide- to a transition-level table.
#'
#' If the SWATH data should be analyzed on transition-level the data needs to be
#' tranformed from peptide-level table to a transition-level table (one row per
#' transition instead of one row per peptide). The columns
#' "aggr_Fragment_Annotation" and "aggr_Peak_Area" are disaggregated into the
#' new columns "Fragmentation" and "Intensity".
#' The following columns are renamed if they exist: FullPeptideName ->
#' PeptideSequence, Charge -> PrecursorCharge, Area -> Intensity, Fragment ->
#' Fragmentation, Sequence -> NakedSequence.
#'
#' @param data A data frame containing SWATH data.
#' @param all.columns Option that all columns are processed. Otherwise only the
#'   columns typically needed for downstream analysis are processed.
#' @return Returns a data frame containing the SWATH data in a transition-level
#'   table.
#' @author Peter Blattmann
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.filtered.decoy <- filter_mscore(data, 0.01)
#'  raw <- disaggregate(data.filtered.decoy)
#'  }
#' @importFrom reshape2 colsplit
#' @export
disaggregate <- function(data, all.columns = FALSE) {

    # sanity test on the number of transitions per precursor
    n.transitions <- lapply(as.character(data$aggr_Fragment_Annotation), 
                            function(x) strsplit(x, ";"))
    n.transitions2 <- unlist(lapply(n.transitions, 
                                    function(x) length(unlist(x))))

    n.transitions3 <- lapply(as.character(data$aggr_Peak_Area), 
                             function(x) strsplit(x, ";"))
    n.transitions4 <- unlist(lapply(n.transitions3, 
                                    function(x) length(unlist(x))))

    if (sum(n.transitions2 != n.transitions4) > 0) {
        stop(paste("The number of transitions annotated and measured do not match in the following transitions:\n",
            paste(unlist(n.transitions[n.transitions2 != n.transitions4]), collapse = ", ")))
    }

    # test if always the same number of transitions per precursor were used
    if (min(n.transitions2) == max(n.transitions2)) {
        message(paste("The library contains", max(n.transitions2), "transitions per precursor.
                  \nThe data table was transformed into a table containing one row per transition."))
    }
    if (min(n.transitions2) != max(n.transitions2)) {
        message(paste("The library contains between", min(n.transitions2), "and",
            max(n.transitions2), "transitions per precursor.\nThe data table was transformed into a table containing one row per transition."))
    }

    data.new <- cbind(data, 
                      colsplit(data$aggr_Fragment_Annotation, ";", paste("Split_FragAnnot_", seq_len(max(n.transitions2)), sep = "")), 
                      colsplit(data$aggr_Peak_Area, ";", paste("Split_PeakArea_", seq_len(max(n.transitions2)), sep = "")), stringsAsFactors = FALSE)

    data.new.m <- reshape2::melt(data.new, 
                                 id.vars = grep("Split_FragAnnot", colnames(data.new),invert = TRUE), 
                                 measure.vars = grep("Split_FragAnnot", colnames(data.new)),
                                 variable.name = "FragAnnot_N", value.name = "Fragment")
    
    data.new.m2 <- reshape2::melt(data.new.m, 
                                  id.vars = grep("Split_PeakArea", colnames(data.new.m), invert = TRUE),
                                  measure.vars = grep("Split_PeakArea", colnames(data.new.m)),
                                  variable.name = "Area_N", value.name = "Area")

    # added because it didn't name the variables in a later trial
    if (sum(colnames(data.new.m2) %in% c("FragAnnot_N")) == 0) {
        l <- length(colnames(data.new.m2))
        colnames(data.new.m2)[l - 3] <- "FragAnnot_N"
        colnames(data.new.m2)[l - 2] <- "Fragment"
        colnames(data.new.m2)[l - 1] <- "Area_N"
        colnames(data.new.m2)[l] <- "Area"
    }


    data.new.m3 <- data.new.m2[gsub("Split_FragAnnot_", "", 
                                    data.new.m2[, "FragAnnot_N"]) == gsub("Split_PeakArea_", "", data.new.m2[, "Area_N"]), ]

    if (!isTRUE(all.columns)) {
        cols <- colnames(data.new.m3)[colnames(data.new.m3) %in% c("ProteinName",
            "FullPeptideName", "PeptideSequence", "Sequence", "Charge", 
            "PrecursorCharge", "Fragment", "FragmentIon", "Area", "Condition", 
            "BioReplicate", "Run", "RT")]
        data.new.merged <- data.new.m3[, cols]
    }
    if (isTRUE(all.columns)) {
        data.new.merged <- data.new.m3
    }

    colnames(data.new.merged) <- gsub("FullPeptideName", "PeptideSequence", 
                                      colnames(data.new.merged))
    colnames(data.new.merged) <- gsub("^Charge$", "PrecursorCharge", 
                                      colnames(data.new.merged))

    colnames(data.new.merged) <- gsub("Area", "Intensity", 
                                      colnames(data.new.merged))
    colnames(data.new.merged) <- gsub("Fragment", "FragmentIon", 
                                      colnames(data.new.merged))

    if ("Sequence" %in% cols) {
        colnames(data.new.merged) <- gsub("^Sequence$", "NakedSequence", 
                                          colnames(data.new.merged))
    }
    if (sum(is.na(data.new.merged$Intensity)) > 0) {
        .ids <- !is.na(data.new.merged$Intensity)
        message(paste((length(data.new.merged$Intensity) - sum(.ids)), "row(s) was/were removed because they did not contain data due to different number of transitions per precursor"))
        data.new.merged <- data.new.merged[.ids, ]
    }

    return(data.new.merged)
}
