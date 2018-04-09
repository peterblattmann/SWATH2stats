sample_annotation <- function(data, sample_annotation, data_type="OpenSWATH",
                              annotation_file_column="Filename",
                              data_file_column="filename",
                              condition_column="Condition",
                              replicate_column="BioReplicate",
                              run_column="Run",
                              change_run_id=TRUE, verbose=FALSE) {
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
      warning(paste0("The missing data samples from the annotation are: ", missing_samples_from_data, "."))
    }
    if (missing_samples_from_annot > 0) {
      warning(paste0("The missing data samples from the data are: ", missing_samples_from_annotation, "."))
    }
  }

  for (i in 1:nrow(sample_annotation)) {
    ## Why not just do that as x %in% y?
    n_found <- grep(pattern=sample_annotation[i, annotation_file_column],
                    x=sample_annotation[, annotation_file_column])
    if (length(n_found) > 1) {
      stop("The values in the column filename are not unique and will lead to erroneous results because the following string matches to multiple different filenames: ", sample_annotation[i,"Filename"])
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
