#' I do not yet understand what is going on in this function.  It is querying
#' the swath2stats data for a few parameters and eventually melting
#' it. Depending on what comes of that, it drops some rows and returns a new,
#' presumably better data structure.
#'
#' @param data SWATH2stats data structure to beat on.
#' @param all.columns  Mess with the set of returned columns?
#' @return a new version of the swath2stats data.
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

  data.new <- cbind(data, colsplit(data[["aggr_fragment_annotation"]], ";",
                                   paste0("Split_FragAnnot_", seq(1:max(n.transitions2)))),
                    colsplit(data[["aggr_peak_area"]], ";",
                             paste0("Split_PeakArea_",
                                    seq(1:max(n.transitions2)))), stringsAsFactors=FALSE)

  data.new.m <- reshape2::melt(data.new, id.vars=grep("Split_FragAnnot",
                                                      colnames(data.new), invert=TRUE),
                               measure.vars=grep("Split_FragAnnot", colnames(data.new)),
                               variable.name="FragAnnot_N", value.name = "Fragment")
  data.new.m2 <- reshape2::melt(data.new.m, id.vars=grep("Split_PeakArea",
                                                         colnames(data.new.m), invert=TRUE),
                                measure.vars=grep("Split_PeakArea", colnames(data.new.m)),
                                variable.name="Area_N", value.name = "Area")

  # added because it didn't name the variables in a later trial
  if (sum(colnames(data.new.m2) %in% c("FragAnnot_N")) == 0) {
    l <- length(colnames(data.new.m2))
    colnames(data.new.m2)[l-3] <- "FragAnnot_N"
    colnames(data.new.m2)[l-2] <- "Fragment"
    colnames(data.new.m2)[l-1] <- "Area_N"
    colnames(data.new.m2)[l] <- "Area"
  }


  data.new.m3 <- data.new.m2[gsub("Split_FragAnnot_", "",
                                  data.new.m2[,"FragAnnot_N"]) == gsub("Split_PeakArea_", "",
                                                                       data.new.m2[,"Area_N"]), ]

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
            "row(s) was/were removed because they did not contain data due to different number of transitions per precursor")
    data.new.merged <- data.new.merged[.ids, ]
  }

  return(data.new.merged)
}
