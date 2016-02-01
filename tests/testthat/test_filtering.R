context("data processing")

test_that("filtering", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data(Study_design, package="SWATH2stats")

  data <- reduce_OpenSWATH_output(data)
  data <- sample_annotation(data, Study_design)

  ## assess_decoy_rate
  assess_decoy_rate(data)
  expect_that(assess_decoy_rate(data), shows_message("Number of non-decoy peptides: 273"))
  expect_that(assess_decoy_rate(data), shows_message("Number of decoy peptides: 11"))
  expect_that(assess_decoy_rate(data), shows_message("Decoy rate: 0.0403"))

  ## filter_mscore
  data.filtered.mscore <- filter_mscore(data, 0.01)
  data.filtered.mscore <- filter_mscore(data, 0.001)
  data.filtered.mscore <- filter_mscore(data, 0.0001)

  expect_that(filter_mscore(data, 0.01), shows_message("Dimension difference: 931, 0"))
  expect_that(filter_mscore(data, 0.001), shows_message("Dimension difference: 1068, 0"))
  expect_that(filter_mscore(data, 0.0001), shows_message("Dimension difference: 1166, 0"))

  # filter mscore_requant
  data.filtered.mscore <- filter_mscore_freqobs(data, 0.01, 0.8)

  expect_that(filter_mscore_freqobs(data, 0.01, 0.8), shows_message("Treshold, peptides need to have been quantified in more conditions than: 4.8"))
  expect_that(filter_mscore_freqobs(data, 0.01, 0.8), shows_message("Fraction of peptides selected: 0.42"))
  expect_that(filter_mscore_freqobs(data, 0.01, 0.8), shows_message("Dimension difference: 1323, 0"))

  # filter mscore_conditions
  data.filtered.mscore <- filter_mscore_condition(data, 0.01, 3)
  expect_that(filter_mscore_condition(data, 0.01, 3), shows_message("Fraction of peptides selected: 0.47"))
  expect_that(filter_mscore_condition(data, 0.01, 3), shows_message("Dimension difference: 1209, 0"))

  # filter on_fdr
  expect_that(filter_mscore_fdr(data, FFT = 0.7, overall_protein_fdr_target=0.1, upper_overall_peptide_fdr_limit=0.05),
  shows_message("been removed from the returned data"))
})

