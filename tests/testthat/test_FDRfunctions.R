context("FDR assessment on annotated data")

test_that("FDR assessment on annotated data", {
  data(OpenSWATH_data, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")
  data.FDR<-sample_annotation(OpenSWATH_data, Study_design)

  # assess_fdr_overall output test
  expect_message(assess_fdr_overall(data.FDR), ".csv")
  expect_message(assess_fdr_overall(data.FDR), ".pdf")
  expect_that(length(assess_fdr_overall(data.FDR)), equals(0))
  expect_that(length(assess_fdr_overall(data.FDR, output = "Rconsole")), equals(13))
  
  # assess_fdr_byrun output test
  expect_message(assess_fdr_byrun(data.FDR), ".csv")
  expect_message(assess_fdr_byrun(data.FDR), ".pdf")
  expect_message(assess_fdr_byrun(data.FDR), "assay")
  expect_message(assess_fdr_byrun(data.FDR), "peptide")
  expect_message(assess_fdr_byrun(data.FDR), "protein")

  # assess_fdr_byrun return test
  expect_that(round(mean(assess_fdr_byrun(data.FDR, output = "Rconsole")[1,,]),digits = 4), equals(229.3333))

  # mscore4(assay,pep,prot)fdr return tests
  expect_that(signif(mscore4assayfdr(data.FDR, fdr_target = 0.001), digits=5), equals(0.00017783))
  expect_that(signif(mscore4pepfdr(data.FDR, fdr_target = 0.002), digits=5), equals(0.00017783))
  expect_that(signif(mscore4protfdr(data.FDR, fdr_target = 0.05), digits=5), equals(0.00017783))

})
