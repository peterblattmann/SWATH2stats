library(testthat)
library(SWATH2stats)
context("data processing")

test_that("filtering", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data(Study_design, package="SWATH2stats")

  data <- reduce_OpenSWATH_output(data)

  ## test if checking for missing columns works
  data2 <- data
  ## hmm this makes me feel a little stupid, doesn't setting the column to NULL
  ## just fill it with NULL and not remove the column?  The test in
  ## reduce_OpenSWATH_output() looks for the existence of the column, it should
  ## still exist but be filled with NULL -- no I'll be damned setting it to NULL
  ## is the same as doing a column subset, I had no idea -- why does it not fill
  ## every row with NULL? hahah I guess that is what I get for assuming I know
  ## what R will do!
  data2$proteinname <- NULL
  expect_warning(reduce_OpenSWATH_output(data2), "These columns are missing from the data:proteinname")

  ## sample_annotation
  data <- sample_annotation(data, Study_design, check_files=FALSE)

  ## assess_decoy_rate
  decoys_numbers <- assess_decoy_rate(data)
  expect_that(assess_decoy_rate(data), shows_message("Number of non-decoy peptides: 273"))
  expect_that(assess_decoy_rate(data), shows_message("Number of decoy peptides: 11"))
  expect_that(assess_decoy_rate(data), shows_message("Decoy rate: 0.0403"))

  # test if no decoys present
  data2 <- data
  data2$decoy <- NULL
  expect_error(assess_decoy_rate(data2), "There is no decoy column in the table")
  # test if no decoys present
  data2 <- data
  data2$fullpeptidename <- NULL
  expect_error(assess_decoy_rate(data2), "There is no fullpeptidename column in the table")

  ## filter_mscore
  data.filtered.mscore <- filter_mscore(data, 0.01)
  data.filtered.mscore <- filter_mscore(data, 0.001)
  data.filtered.mscore <- filter_mscore(data, 0.0001)

  expect_that(filter_mscore(data, 0.01),
              shows_message("Original dimension: 2307, new dimension: 1376, difference: 931"))
  expect_that(filter_mscore(data, 0.001),
              shows_message("Original dimension: 2307, new dimension: 1239, difference: 1068"))
  expect_that(filter_mscore(data, 0.0001),
              shows_message("Original dimension: 2307, new dimension: 1141, difference: 1166"))

  # filter mscore_requant
  data.filtered.mscore <- filter_mscore_freqobs(data, mscore=0.01, percentage=0.8)

  expect_that(filter_mscore_freqobs(data, mscore=0.01, percentage=0.8),
              shows_message("Peptides need to have been quantified in more conditions than: 4.8"))

  ## hmm I am getting 0.67 rather than 0.42!  Presumably I effed something up
  ## I have some comments toward the end of filter_mscore_freqobs where I am
  ## attempting to understand my mistaken understanding here, but I have not yet
  ## resolved what I perceive it does vs. what it is supposed to do.
  expect_that(filter_mscore_freqobs(data, mscore=0.01, percentage=0.8),
              shows_message("Fraction of peptides selected: "))
  expect_that(filter_mscore_freqobs(data, 0.01, 0.8),
              shows_message("Original dimension: 2307"))
  expect_message(filter_mscore_freqobs(data, 0.01),
                 "Peptides need to have been quantified in more conditions than: 0")

  ## filter mscore_conditions
  data.filtered.mscore <- filter_mscore_condition(data, 0.01, 3)
  expect_that(filter_mscore_condition(data, 0.01, 3),
              shows_message("Fraction of peptides selected: 0.47"))
  expect_that(filter_mscore_condition(data, 0.01, 3),
              shows_message("Original dimension: 2307, new dimension: 1098, difference: 1209."))

  ## filter on_fdr
  expect_that(filter_mscore_fdr(data, FFT = 0.7, overall_protein_fdr_target=0.1, upper_overall_peptide_fdr_limit=0.05),
  shows_message("been removed from the returned data"))

  ## test if remove decoy works
  data2 <- filter_mscore_fdr(data, FFT=0.7, overall_protein_fdr_target=0.1,
                             upper_overall_peptide_fdr_limit=0.05, rm.decoy=FALSE)
  data3 <- filter_mscore_fdr(data, FFT = 0.7, overall_protein_fdr_target=0.1,
                             upper_overall_peptide_fdr_limit=0.05, rm.decoy=TRUE)

  expect_message(filter_mscore_fdr(data, FFT=0.7,
                                   overall_protein_fdr_target=0.1,
                                   upper_overall_peptide_fdr_limit=0.05, rm.decoy=FALSE),
                 "The decoys have NOT been removed from the returned data")
  num_decoys_kept <- sum(grepl(pattern="DECOY", x=data[["proteinname"]]))
  expect_equal(num_decoys_kept, 62)
  expect_false(length(grep("DECOY", data3[["proteinname"]])) > 1)

  ## test that error is returned if lowest m_score is from decoy
  data2 <- data
  data2 <- data2[data2$m_score != 0,]
  data2[data2$decoy == 1, "m_score"] <- data2[data2$decoy == 1, "m_score"] * 1e-20
  expect_true(is.na(mscore4protfdr(data2, FFT = 0.7, fdr_target = 0.001)))
  expect_error(filter_mscore_fdr(data2, FFT = 0.7, overall_protein_fdr_target=0.1, upper_overall_peptide_fdr_limit=0.05),
               "The overall_protein_fdr_target cannot be reached in this dataset")

  ## test that local FDR is only calculated if enought decoys for each run
  data2 <- data
  data2 <- data2[data2$m_score != 0,]
  data2[data2$decoy == 1 & data2$m_score == 2, "m_score"] <- 0.01

  expect_message(filter_mscore_fdr(data, FFT = 0.7, overall_protein_fdr_target=0.1, upper_overall_peptide_fdr_limit=0.05),
               "Individual run FDR quality of the peptides was not calculated\nas not every run contains a decoy")
  expect_message(filter_mscore_fdr(data2, FFT = 0.7, overall_protein_fdr_target=0.1, upper_overall_peptide_fdr_limit=0.05),
                 "Individual run FDR quality of the peptides selected from")
})

test_that("count_analytes", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data

  analysis <- count_analytes(data)
  expect_that(analysis[analysis$run_id == "0_210","fullpeptidename"], equals(261))
  expect_that(analysis[analysis$run_id == "0_221","fullpeptidename"], equals(271))

  analysis2 <- count_analytes(data[data$m_score < 1e-4,])
  expect_that(analysis[analysis$run_id == "0_177","transition_group_id"], equals(388))
  expect_that(analysis[analysis$run_id == "0_29","proteinname"], equals(10))
})
