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
#' @param data  A data frame containing SWATH data.
#' @param sample_annotation A data frame containing the columns: Filename,
#'   Condition, BioReplicate, Run. The values contained in the column filename
#'   have to be present in the filename of the SWATH data.
#' @param data_type Option to specify the format of the table, if the column
#'   names from an OpenSWATH output or MSstats table are used.
#' @param annotation_file_column  Name of the column containing the output file
#'   from the experiment metadata for  each sample.  In my own sample sheet, I keep
#'   columns for the mzXML files, tric outputs, raw files, and osw outputs from
#'   OpenSwathWorkFlow; and I cannot be relied upon to remember which is which,
#'   ergo this option.
#' @param data_file_column Option to specify the column name where the injection
#'   file is specified. Default is set to "filename".
#' @param condition_column  Which column annotates the experimental condition in
#'   the swath data?
#' @param replicate_column  Which column annotates the replicate in the data?
#' @param run_column  Which column annotates the separate runs?
#' @param change_run_id  Option to choose if the run\_id column shall be
#'   reassigned to a unique value combining the values of Condition,
#'   BioReplicate and Run. (Option only possible if data is of format
#'   "OpenSWATH")
#' @param verbose  Option to turn on reporting on which filename it is working
#'   on.
#' @return Returns a dataframe with each row annotated for the study design
#' @author Peter Blattman
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- SWATH2stats::sample_annotation(OpenSWATH_data, Study_design, verbose=TRUE)
#'  summary(data)
#' @export
sample_annotation <- function(data, sample_annotation, data_type="OpenSWATH",
                              annotation_file_column="filename",
                              check_files=TRUE,
                              data_file_column="filename",
                              condition_column="condition",
                              replicate_column="bioreplicate",
                              fullpeptidename_column=c("fullpeptidename", "fullunimodpeptidename"),
                              run_id=NULL,
                              run_column="run",
                              change_run_id=TRUE,
                              verbose=FALSE) {
  #### annotate sample
  ### needs a txt file with the columns Filename, Condition, BioReplicate, Run. In Filename a unique string contained in File
  ### must be contained.
  colnames(sample_annotation) <- tolower(colnames(sample_annotation))
  colnames(data) <- tolower(colnames(data))
  data_type <- tolower(data_type)
  annotation_file_column <- tolower(annotation_file_column)
  data_file_column <- tolower(data_file_column)
  condition_column <- tolower(condition_column)
  replicate_column <- tolower(replicate_column)
  run_column <- tolower(run_column)
  fullpeptidename_column <- tolower(fullpeptidename_column)

  if (! class(data)[1] == "data.frame" | class(data)[1] == "data.table") {
    stop("Input data is not a data.frame/data.table")
  }
  if (! annotation_file_column %in% colnames(sample_annotation)) {
    stop(paste0("The file column: '", annotation_file_column, "' is missing from the annotations."))
  }
  if (! data_file_column %in% colnames(data)) {
    stop(paste0("The data file column: '", data_file_column, "' is missing from the data file."))
  }

  annotation_files <- levels(as.factor(sample_annotation[[annotation_file_column]]))
  data_files <- levels(as.factor(data[[data_file_column]]))
  ## I want to make sure all these functions work for the data provided with the package.
  ## Unfortunately, the files listed in the data look like:
  ## '/scratch/9148912somethingsomething/filename.mzXML.gz' while the files in
  ## the annotations look like 'filename'
  ## Therefore my test that the two sets are equivalent is a priori doomed to
  ## failure.
  annotation_files <- basename(annotation_files)
  data_files <- basename(data_files)
  removal_regexes <- c("\\.gz$", "\\.mzXML$", "\\.mzML$", "\\.xls", "\\.xlsx")
  for (removal in removal_regexes) {
    annotation_files <- gsub(pattern=removal, replacement="",
                             x=annotation_files)
    data_files <- gsub(pattern=removal, replacement="",
                       x=data_files)
  }
  annotation_files <- sort(annotation_files)
  data_files <- sort(data_files)

  if (isTRUE(check_files)) {
    if (isTRUE(all.equal(annotation_files, data_files))) {
      message("Found the same mzXML files in the annotations and data.")
    } else {
      warning("The files listed in the sample annotation and data are not equivalent.")
      missing_samples_from_data_idx <- ! annotation_files %in% data_files
      missing_samples_from_data <- annotation_files[missing_samples_from_data_idx]
      missing_samples_from_annot_idx <- ! data_files %in% annotation_files
      missing_samples_from_annot <- data_files[missing_samples_from_annot_idx]
      if (missing_samples_from_data > 0) {
        warning(paste0("The missing data samples from the annotation are: ",
                       toString(missing_samples_from_data), "."))
      }
      if (missing_samples_from_annot > 0) {
        warning(paste0("The missing data samples from the data are: ",
                       toString(missing_samples_from_annot), "."))
      }
    }
  } else {
    message("Not checking that the files are identical between the annotation and data.")
  }

  for (i in seq_len(nrow(sample_annotation))) {
    ## Why not just do that as x %in% y?
    n_found <- grep(pattern=sample_annotation[i, annotation_file_column],
                    x=sample_annotation[, annotation_file_column])
    if (length(n_found) > 1) {
      stop("The values in the column filename are not unique and will lead to erroneous results because the following string matches to multiple different filenames: ", sample_annotation[i, "Filename"])
    }
  }

  if (! data_type %in% c("openswath", "msstats")) {
    ## Only accept openswath or msstats inputs
    stop("Type of data not recognized")
  } else {
    for (i in levels(factor(sample_annotation[[annotation_file_column]]))) {
      if (verbose) {
        message(i)
      }
      coord <- grep(i, data[, data_file_column], fixed=TRUE)
      if (length(coord) == 0) {
        warning(paste0("No measurement value found for sample '", i, "' in the data file."))
      }

      data_subset <- sample_annotation[which(i == sample_annotation[[annotation_file_column]]), ]
      data[coord, "condition"] <- data_subset[, condition_column]
      data[coord, "bioreplicate"] <- data_subset[, replicate_column]
      data[coord, "run"] <- data_subset[, run_column]
    }

    # select column names
    if (data_type == "openswath") {
      full_peptide <- fullpeptidename_column %in% colnames(data)
      if (sum(full_peptide) == 0) {
        stop("Unable to find a column containing the full peptide name.")
      } else if (sum(full_peptide) == 1) {
        full_peptide <- fullpeptidename_column[full_peptide]
      } else {
        chosen_full_peptide <- fullpeptidename_column[full_peptide][1]
        warning("There are multiple possibilities for the full peptidename, ",
                chosen_full_peptide, " was chosen.")
        full_peptide <- chosen_full_peptide
      }
      add_colnames <- colnames(data)[!(colnames(data) %in%
                                       c("proteinname", full_peptide, "charge",
                                         "aggr_fragment_annotation", "aggr_peak_area",
                                         "condition", "bioreplicate", "run"))]

      sel_colnames <- c("proteinname", full_peptide, "charge",
                        "aggr_fragment_annotation", "aggr_peak_area",
                        "condition", "bioreplicate", "run", add_colnames)

      missing_columns <- ! sel_colnames %in% colnames(data)
      ## Oh come _on_ tric change the capitalization of some columns!
      ## I _really_ resented MSstats' choice to do toupper() on all the data
      ## columns, but it begins to make sense...
      colnames(data) <- tolower(colnames(data))
      sel_colnames <- tolower(sel_colnames)
      if (sum(missing_columns) > 0) {
        message("The following columns were missing from the data:")
        message(toString(sel_colnames[missing_columns]))
      }
      data <- data[, sel_colnames]
      colnames(data) <- c("proteinname", "fullpeptidename", "charge",
                          "aggr_fragment_annotation", "aggr_peak_area",
                          "condition", "bioreplicate", "run", add_colnames)
    }

    if (data_type == "msstats") {
      add_colnames <- colnames(data)[!(colnames(data) %in%
                                       c("proteinname", "peptidesequence", "precursorcharge",
                                         "fragmention", "productcharge", "isotopelabeltype",
                                         "condition", "bioreplicate", "run", "intensity"))]

      sel_colnames <- c("proteinname", "peptidesequence", "precursorcharge", "fragmention",
                        "productcharge", "isotopelabeltype", "condition", "bioreplicate",
                        "run", "intensity", add_colnames)
      data <- data[, sel_colnames]
    }

    if (isTRUE(change_run_id)) {
      ##if (is.null(run_column)) {
        data[["run_id"]] <- paste(data[["condition"]],
                                  data[["bioreplicate"]],
                                  data[["run"]], sep="_")
      ##} else {
      ##  data[["run_id"]] <- data[[run_column]]
      ##} ## Use a specific column for the run id?
    }  ## change the run id?
  }  ### OpenSwath data
  if (isTRUE(verbose)) {
    message(nrow(sample_annotation), " samples were read from the annotations.")
    message(nrow(data), " transitions were read from the data and merged with the annotations.")
  }
  return(data)
}
