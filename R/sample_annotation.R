#' Given dataframes of TRIC processed data and sample annotations, mash them
#' together into something appropriate for downstream analyses.
#'
#' This performs some quick sanity checks on the data and annotations and
#' creates the 'Condition', 'BioReplicate', and 'Run' columns along with other
#' columns expected by MSstats/OpenSWATH.
#'
#' @param data  Data frame following processing by feature_alignment.py.
#' @param data_type  Either OpenSWATH or MSstats, this defines the columns which
#'   get created in the output data frame.
#' @param annotation_file_column  Name of the column containing the output file
#'   from the experiment metadata for  each sample.  In my own sample sheet, I keep
#'   columns for the mzXML files, tric outputs, raw files, and osw outputs from
#'   OpenSwathWorkFlow; and I cannot be relied upon to remember which is which,
#'   ergo this option.
#' @param data_file_column  Column in the feature_alignment.py containing the
#'   mzXML file for each detected transition.
#' @param condition_column  Which column annotates the experimental condition in
#'   the swath data?
#' @param replicate_column  Which column annotates the replicate in the data?
#' @param run_column  Which column annotates the separate runs?
#' @param change_run_id  Concatenate the condition, replicate, and run columns
#'   into a new run_id column?
#' @param verbose  Be verbose?
#' @return Dataframe with each row annotated from the sample metadata.
#' @author Peter Blattman
#' @examples
#'  \dontrun{
#'    data("OpenSWATH_data", package="SWATH2stats")
#'    data("Study_design", package="SWATH2stats")
#'    data <- sample_annotation(OpenSWATH_data, Study_design)
#' }
#' @export
sample_annotation <- function(data, sample_annotation, data_type="OpenSWATH",
                              annotation_file_column="Filename",
                              data_file_column="filename",
                              condition_column="Condition",
                              replicate_column="BioReplicate",
                              run_column="Run", change_run_id=TRUE, verbose=FALSE) {
  #### annotate sample
  ### needs a txt file with the columns Filename, Condition, BioReplicate, Run. In Filename a unique string contained in File
  ### must be contained.
  if (! class(data) == "data.frame" | class(data) == "data.table") {
    stop("Input data is not a data.frame/data.table")
  }
  if (! annotation_file_column %in% colnames(sample_annotation)) {
    stop(paste0("The file column: ", annotation_file_column, " is missing from the annotations."))
  }
  if (! data_file_column %in% colnames(data)) {
    stop(paste0("The file column: ", data_file_column, " is missing from the data file."))
  }

  annotation_files <- sort(levels(factor(sample_annotation[[annotation_file_column]])))
  data_files <- sort(levels(factor(data[[data_file_column]])))
  if (! all.equal(annotation_files, data_files)) {
    warning("The number of sample annotation condition and filenames in data are not balanced.")
    missing_samples_from_data_idx <- ! annotation_files %in% data_files
    missing_samples_from_data <- annotation_files[missing_samples_from_data_idx]
    missing_samples_from_annot_idx <- ! data_files %in% annotation_files
    missing_samples_from_annot <- data_files[missing_samples_from_annot_idx]
    if (missing_samples_from_data > 0) {
      warning(paste0("The missing data samples from the annotation are: ",
                     missing_samples_from_data, "."))
    }
    if (missing_samples_from_annot > 0) {
      warning(paste0("The missing data samples from the data are: ",
                     missing_samples_from_annotation, "."))
    }
  }

  for (i in 1:nrow(sample_annotation)) {
    ## Why not just do that as x %in% y?
    n_found <- grep(pattern=sample_annotation[i, annotation_file_column],
                    x=sample_annotation[, annotation_file_column])
    if (length(n_found) > 1) {
      stop("The values in the column filename are not unique and will lead to erroneous results because the following string matches to multiple different filenames: ", sample_annotation[i, "Filename"])
    }
  }

  if (data_type %in% c("OpenSWATH", "MSstats")) {
    for (i in levels(factor(sample_annotation[[annotation_file_column]]))) {
      if (verbose) {
        message(i)
      }
      coord <- grep(i, data[, data_file_column], fixed=TRUE)
      if (length(coord) == 0) {
        warning(paste0("No measurement value found for this sample in the data file: ", i, "."))
      }

      data_subset <- sample_annotation[which(i == sample_annotation[[annotation_file_column]]), ]
      data[coord, "Condition"] <- data_subset[, condition_column]
      data[coord, "BioReplicate"] <- data_subset[, replicate_column]
      data[coord, "Run"] <- data_subset[, run_column]
    }

    # select column names
    if (data_type == "OpenSWATH") {
      add_colnames <- colnames(data)[!(colnames(data) %in%
                                       c("ProteinName", "FullPeptideName", "Charge",
                                         "aggr_Fragment_Annotation", "aggr_Peak_Area",
                                         "Condition", "BioReplicate", "Run"))]

      sel_colnames <- c("ProteinName", "FullPeptideName", "Charge", "aggr_Fragment_Annotation",
                        "aggr_Peak_Area", "Condition", "BioReplicate", "Run", add_colnames)
    }
    if (data_type == "MSstats") {
      add_colnames <- colnames(data)[!(colnames(data) %in%
                                       c("ProteinName", "PeptideSequence", "PrecursorCharge",
                                         "FragmentIon", "ProductCharge", "IsotopeLabelType",
                                         "Condition", "BioReplicate", "Run", "Intensity"))]

      sel_colnames <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                        "ProductCharge", "IsotopeLabelType", "Condition", "BioReplicate",
                        "Run", "Intensity", add_colnames)
    }
    data <- data[, sel_colnames]

    if (change_run_id) {
      data[["run_id"]] <- paste(data[["Condition"]],
                                data[["BioReplicate"]],
                                data[["Run"]], sep="_")
    }
    return(data)
  } else {
    message("Type of data not recognized")
  }
}
