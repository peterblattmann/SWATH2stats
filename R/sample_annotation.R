#' Annotate the SWATH data with the sample information
#'
#' For statistical analysis and filtering the measurements need to be annotated
#' with Filename, Condition, BioReplicate, and Run. This functions takes this
#' information from a txt file containing this meta-data.
#'
#' Given dataframes of TRIC processed data and sample annotations, mash them
#' together into something appropriate for downstream analyses.
#'
#' This performs some quick sanity checks on the data and annotations and
#' creates the 'Condition', 'BioReplicate', and 'Run' columns along with other
#' columns expected by MSstats/OpenSWATH.
#'
#' @param data A data frame containing SWATH data.
#' @param sample_annotation A data frame containing the columns: Filename,
#'   Condition, BioReplicate, Run. The values contained in the column filename
#'   have to be present in the filename of the SWATH data.
#' @param data_type Option to specify the format of the table, if the column
#'   names from an OpenSWATH output or MSstats table are used.
#' @param annotation_file_column Name of the column containing the output file
#'   from the experiment metadata for  each sample.  In my own sample sheet, I keep
#'   columns for the mzXML files, tric outputs, raw files, and osw outputs from
#'   OpenSwathWorkFlow; and I cannot be relied upon to remember which is which,
#'   ergo this option.
#' @param check_files Boolean checking if one wishes to ensure that the files
#'   listed in the annotation data are equivalent to the files in the actual
#'   data.
#' @param data_file_column Option to specify the column name where the injection
#'   file is specified. Default is set to "filename".
#' @param condition_column Which column annotates the experimental condition in
#'   the swath data?
#' @param replicate_column Which column annotates the replicate in the data?
#' @param fullpeptidename_column Character list of possible column names used
#'   to define the full peptide name.
#' @param run_column Which column annotates the separate runs?
#' @param change_run_id Option to choose if the run\_id column shall be
#'   reassigned to a unique value combining the values of Condition,
#'   BioReplicate and Run. (Option only possible if data is of format
#'   "OpenSWATH")
#' @param verbose Option to turn on reporting on which filename it is working
#'   on.
#' @return Returns a dataframe with each row annotated for the study design
#' @author Peter Blattman
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- SWATH2stats::sample_annotation(OpenSWATH_data, Study_design, verbose=TRUE)
#'  summary(data)
#' @export
sample_annotation <- function(data, sample.annotation, data.type = "OpenSWATH",
                              column.file = "filename",
                              change.run.id = TRUE, verbose = FALSE) {
    ## annotate sample needs a txt file with the columns Filename, Condition,
    ## BioReplicate, Run. In Filename a unique string contained in File must be
    ## contained.
    if (!isTRUE(class(data) == "data.frame")) {
        stop("Input data is not a data.frame")
    }

    if (!(column.file %in% colnames(data))) {
        stop("Column for filename is not present in data file")
    }
    if (nlevels(factor(sample.annotation$Filename)) != nlevels(factor(data[, column.file]))) {
        warning("The number of sample annotation condition and filenames in data are not balanced.",
            "\n", "Different filenames in sample annotation file: ", nlevels(factor(sample.annotation$Filename)),
            "\n", "Different filenames in data file: ", nlevels(factor(data[, column.file])))
    }

    for (i in seq_len(nrow(sample.annotation))) {
        n.found <- grep(sample.annotation[i, "Filename"], sample.annotation[, "Filename"])
        if (length(n.found) > 1) {
            stop("The values in the column filename are not unique and will lead to erroneous results because the following string matches to multiple different filenames: ",
                sample.annotation[i, "Filename"])
        }
    }

    if (data.type %in% c("OpenSWATH", "MSstats")) {

        for (i in levels(factor(sample.annotation$Filename))) {
            if (verbose) {
                message(i)
            }

            coord <- grep(i, data[, column.file], fixed = TRUE)

            if (length(coord) == 0) {
                warning("No measurement value found for this sample in the data file: ",
                  print(i))
            }

            data.subset <- sample.annotation[which(i == sample.annotation$Filename),
                ]
            data[coord, "Condition"] <- data.subset[, "Condition"]
            data[coord, "BioReplicate"] <- data.subset[, "BioReplicate"]
            data[coord, "Run"] <- data.subset[, "Run"]
        }

        # select column names
        if (data.type == "OpenSWATH") {
            add.colnames <- colnames(data)[!(colnames(data) %in% c("ProteinName",
                "FullPeptideName", "Charge", "aggr_Fragment_Annotation", "aggr_Peak_Area",
                "Condition", "BioReplicate", "Run"))]

            sel.colnames <- c("ProteinName", "FullPeptideName", "Charge", "aggr_Fragment_Annotation",
                "aggr_Peak_Area", "Condition", "BioReplicate", "Run", add.colnames)
        }
        if (data.type == "MSstats") {
            add.colnames <- colnames(data)[!(colnames(data) %in% c("ProteinName",
                "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge",
                "IsotopeLabelType", "Condition", "BioReplicate", "Run", "Intensity"))]

            sel.colnames <- c("ProteinName", "PeptideSequence", "PrecursorCharge",
                "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition",
                "BioReplicate", "Run", "Intensity", add.colnames)
        }
        data <- data[, sel.colnames]

        if (change.run.id) {
            data$run_id <- paste(data$Condition, data$BioReplicate, data$Run, sep = "_")
        }

        return(data)
    } else (print("Type of data not recognized"))
}
