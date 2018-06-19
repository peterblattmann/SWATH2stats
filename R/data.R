#' Sample data from the OpenMS, pyprophet workflow.
#'
#' I am not sure where this data came from, I need to look a little more closely.
#'
#' @usage data(OpenSWATH_data)
#' @format a data frame of 2369 rows and 58 columns.
#' @keywords datasets
#' @examples
#'  data(OpenSWATH_data)
#'  dim(OpenSWATH_data)
"OpenSWATH_data"

#' 6 rows of sample MSstats data.
#'
#' This provides 6 intensities of putative MSstats data for 3 replicates of HeLa cells.
#'
#' @docType data
#'
#' @usage data(MSstats_data)
#'
#' @format a data frame of 6 rows and 11 columns.
#'
#' @keywords datasets
#'
#' @examples
#'  data(MSstats_data)
#'  dim(MSstats_data)
"MSstats_data"

#' A set of Streptococcus pyogenes intensities from 4 samples.
#'
#' To my eyes, this appears to be the result of OpenMS->pyprophet with the
#' exception that the alignment files look to be gzipped Excel files.  I bet
#' they are actually tab separated outputs from feature_alignment.py or
#' pyprophet.
#'
#' @docType data
#'
#' @usage data(Spyogenes)
#'
#' @format a data frame of 32,272 rows and 13 columns.
#'
#' @keywords datasets
#'
#' @examples
#'  data(Spyogenes)
#'  data <- Spyogenes
#'  dim(data)
"Spyogenes"

#' A set of sample annotations which describe an experiment 6 of HeLa control
#' and treated samples (3 of each).
#'
#' This provides Filename, Condition, BioReplicate, and Run data for a small
#' example dataset.  I assumed at first that these correspond to the data in
#' MSstats_data, but that cannot be, as they have different Run IDs.
#'
#' @docType data
#'
#' @usage data(Study_design)
#'
#' @format a data frame of 6 rows and 4 columns.
#'
#' @keywords datasets
#'
#' @examples
#'  data(Study_design)
#'  dim(Study_design)
"Study_design"
