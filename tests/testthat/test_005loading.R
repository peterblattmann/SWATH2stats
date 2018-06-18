context("load data and annotate")

test_that("load data and annotate", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data(Study_design, package="SWATH2stats")

  # filter_openSWATH_output
  expect_that(dim(data)[2], equals(58))
  data <- reduce_OpenSWATH_output(data)
  expect_that(dim(data)[2], equals(14))

  # annotate data
  expect_that(dim(data)[2], equals(14))
  data <- sample_annotation(data, Study_design)

  ## filename.data <- unique(subset(data, Condition == "Hela_Treatment" & BioReplicate == 2 & Run == 4)[,"filename"])
  filename_data_idx <- data[["condition"]] == "Hela_Treatment" &
    data[["bioreplicate"]] == 2 &
    data[["run"]] == 4
  filename_data <- sort(unique(data[filename_data_idx, "filename"]))
  ## filename.design <- subset(Study_design, Condition == "Hela_Treatment" & BioReplicate == 2 & Run == 4)[,"Filename"]
  filename_design_idx <- Study_design[["Condition"]] == "Hela_Treatment" &
    Study_design[["BioReplicate"]] == 2 &
    Study_design[["Run"]] == 4
  filename_design <- sort(unique(Study_design[filename_design_idx, "Filename"]))
  expect_true(length(grep(filename_design, filename_data)) > 0)

  filename_design2 <- subset(Study_design, Condition == "Hela_Treatment" &
                                           BioReplicate == 1 &
                                           Run == 4)[, "Filename"]
  expect_error(length(grep(filename_design2, filename_data)) > 0)

  Study_design2 <- Study_design
  Study_design2$Filename <- gsub("(peterb_L[[:digit:]]{6}).*", "\\1", Study_design2$Filename)

  # error that it has different number of file names
  expect_error(sample_annotation(data, Study_design2))

  # error that the strings in Filename are not unique
  Study_design3 <- Study_design
  Study_design3$Filename <- gsub("peterb_L150514_001_SW", "peterb_L150514", Study_design3$Filename)

  expect_error(sample_annotation(data, Study_design3))

  data <- sample_annotation(data, Study_design, verbose=TRUE)
  expect_message(sample_annotation(data, Study_design, verbose=TRUE), "peterb_J131223_043")
  expect_error(sample_annotation(data, Study_design, data_file_column="Files"),
               "The data file column: 'files' is missing from the data file.")
})

test_that("load MSstats data and annotate", {
  data(MSstats_data, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")

  MSstats_data2 <- MSstats_data
  MSstats_data2$mscore <- 0.01
  expect_message(transform_MSstats_OpenSWATH(MSstats_data), "No column 'mscore' present")
  expect_message(transform_MSstats_OpenSWATH(MSstats_data2), "Additional columns present in the data: filename, mscore")

  expect_warning(sample_annotation(MSstats_data[, c(1:6,9:11)], Study_design,
                                   data_type="MSstats", data_file_column="FileName"),
                 "No measurement value found for sample 'peterb_J131223_054")

  MSstats_data2$FragmentIon <- NULL
  expect_error(transform_MSstats_OpenSWATH(MSstats_data2), "The data frame doesn't contain all required columns.")
})
