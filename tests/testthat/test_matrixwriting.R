context("Overview matrix writing")

test_that("Overview matrix writing", {
  data(OpenSWATH_data_FDR, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")
  data<-sample_annotation(OpenSWATH_data_FDR, Study_design)
  
  # assess_fdr_overall output test
  expect_message(write_matrix_peptides(data), "folder")
  
  # assess_fdr_overall output test
  expect_message(write_matrix_proteins(data), "folder")
    
})
