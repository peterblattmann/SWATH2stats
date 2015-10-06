context("load data and annotate")

test_that("load data and annotate", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data(Study_design, package="SWATH2stats")
  
  # filter_openSWATH_output
  expect_that(dim(data)[2], equals(58))
  data <- reduce_OpenSWATH_output(data)
  expect_that(dim(data)[2], equals(12))
  
  # annotate data
  expect_that(dim(data)[2], equals(12))
  data <- sample_annotation(data, Study_design)
  
  filename.data <- unique(subset(data, Condition == "Hela_Treatment" & BioReplicate == 2 & Run == 4)[,"align_origfilename"])
  filename.design <- subset(Study_design, Condition == "Hela_Treatment" & BioReplicate == 2 & Run == 4)[,"Filename"]
  filename.design2 <- subset(Study_design, Condition == "Hela_Treatment" & BioReplicate == 1 & Run == 4)[,"Filename"]
  
  expect_true(length(grep(filename.design, filename.data))>0)
  expect_error(length(grep(filename.design2, filename.data))>0)
  
})

