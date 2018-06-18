context("data conversion")

test_that("data conversion", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data(Study_design, package="SWATH2stats")

  data <- sample_annotation(data, Study_design)
  data.filtered.mscore <- filter_mscore_freqobs(data, mscore=0.01, percentage=0.8)

  data2 <- data.table::data.table(data)
  expect_error(sample_annotation(data2, Study_design))

  data.proteotypic <- filter_proteotypic_peptides(data.filtered.mscore)
  expect_that(filter_proteotypic_peptides(data.filtered.mscore),
              shows_message("Number of proteins detected: 11"))

  expect_that(filter_proteotypic_peptides(data.filtered.mscore),
                shows_message("Number of proteotypic peptides detected: 131"))
  expect_that(dim(data.proteotypic), equals(c(972, 62)))

  data.all <- filter_all_peptides(data.filtered.mscore)
  expect_that(filter_proteotypic_peptides(data.filtered.mscore),
              shows_message("Number of proteins detected: 11"))
  expect_that(dim(data.all), equals(c(984, 62)))

  data.max <- filter_on_max_peptides(data.filtered.mscore, n_peptides=5)
  expect_that(dim(data.max), equals(c(300, 62)))
  ## An inteesting(annoying) side-note: formatting strings with 2 return
  ## characters fails, ergo me explicitly typing 1 \n below then the return character...
  expect_that(filter_on_max_peptides(data.filtered.mscore, n_peptides=5),
              shows_message("Before filtering:
  Number of proteins: 10
  Number of peptides: 133\n
Percentage of peptides removed: 69.17%\n
After filtering:
  Number of proteins: 10
  Number of peptides: 41"))

  data.min <- filter_on_min_peptides(data.filtered.mscore, 5)
  expect_that(dim(data.min), equals(c(918, 62)))

  data.python <- convert_python(data.max)
  expect_that(identical(grep("UniMod\\:", data.max[, "aggr_fragment_annotation"]),
                        grep("UniMod\\_", data.python[, "aggr_fragment_annotation"])), is_true())
  expect_that(identical(grep("UniMod\\:", data.max[, "aggr_fragment_annotation"]),
                        grep("UniMod\\_", data.max[, "aggr_fragment_annotation"])), is_false())
  expect_that(identical(grep("UniMod\\:", data.max[, "fullpeptidename"]),
                        grep("UniMod\\_", data.python[, "fullpeptidename"])), is_true())
  expect_that(identical(grep("UniMod\\:", data.max[, "fullpeptidename"]),
                        grep("UniMod\\_", data.max[, "fullpeptidename"])), is_false())

  raw <- disaggregate(data.max)
  expect_true(subset(raw, fragmention == "1069078_FIIDPAAVITGR_2" & run == 6)[, "intensity"] == 188569)
  expect_true(subset(raw, fragmention == "279546_FPSIVGVAR_2" & run == 5)[, "intensity"] == 10805)

  # test error if transitions and data don't match
  data.max.test1 <- data.max
  data.max.test1$aggr_peak_area <- as.character(data.max.test1$aggr_peak_area)
  data.max.test1[data.max.test1$run == 6 &
                 data.max.test1$bioreplicate == 3 &
                 data.max.test1$peptide_charge == "FIIDPAAVITGR 2","aggr_peak_area"] <- "188569.000000;81947.000000;77295.000000;50218.000000"
  expect_error(disaggregate(data.max.test1))
  #expect_that(disaggregate(data.max.test1), throws_error("Error in disaggregate(data.max.test1) : \n  The number of transitions annotated and measured do not match in the following transitions:\n 1069078_FIIDPAAVITGR_2, 1069082_FIIDPAAVITGR_2, 1069080_FIIDPAAVITGR_2, 1069071_FIIDPAAVITGR_2, 1069084_FIIDPAAVITGR_2, 1069076_FIIDPAAVITGR_2\n", fixed=TRUE))
  expect_that(disaggregate(data.max.test1),
              throws_error("The number of transitions annotated and measured do not match in the following transitions:\n 1069078_FIIDPAAVITGR_2, 1069082_FIIDPAAVITGR_2, 1069080_FIIDPAAVITGR_2, 1069071_FIIDPAAVITGR_2, 1069084_FIIDPAAVITGR_2, 1069076_FIIDPAAVITGR_2", fixed=FALSE))

  # test if spectral library does not contain same number of transitions
  data.max.test2 <- data.max.test1
  data.max.test2$aggr_fragment_annotation <- as.character(data.max.test2$aggr_fragment_annotation)
  data.max.test2[data.max.test1$run == 6 &
                 data.max.test1$bioreplicate == 3 &
                 data.max.test1$peptide_charge == "FIIDPAAVITGR 2",
                 "aggr_fragment_annotation"] <- "1069078_FIIDPAAVITGR_2;1069082_FIIDPAAVITGR_2;1069080_FIIDPAAVITGR_2;1069071_FIIDPAAVITGR_2"
  expect_that(disaggregate(data.max.test2),
              shows_message("The library contains between 4 and 6 transitions per precursor.", fixed=TRUE))
  expect_that(dim(disaggregate(data.max.test2)), equals(c(1798, 10)))

  data.MSstats <- convert_MSstats(raw)
  expect_that(dim(data.MSstats), equals(c(1800,10)))

  data.mapDIA <- convert_mapDIA(raw, RT = TRUE)
  expect_that(dim(data.mapDIA), equals(c(300,10)))

  data.aLFQ <- convert_aLFQ(raw)
  expect_that(dim(data.aLFQ), equals(c(1800,8)))

  expect_warning(convert_aLFQ(raw[c(1:6),]),
                 "The aLFQ package should only be used with transition-level data.")

  # test if warning is displayed when there are several values for a data point
  raw2 <- raw
  raw2$intensity <- raw$intensity + 20
  raw2 <- rbind(raw[raw$fragmention == "1069078_FIIDPAAVITGR_2" &
                    raw$condition %in% c("Hela_Treatment"),], raw2)
  expect_warning(convert_mapDIA(raw2, RT=TRUE), "Data contains several intensity values per condition")

  # test if warning is displayed when column is missing
  raw2 <- raw
  raw2$rt <- NULL
  expect_warning(convert_mapDIA(raw2, RT=TRUE),
                 "One or several columns required by mDIA were not in the data and filled with NAs")

  # test if extraction of Precursor charge works
  raw2 <- raw
  raw2$precursorcharge <- NULL
  expect_message(convert_MSstats(raw2),
                 "Missing columns: precursorcharge")
  data.MSstats <- convert_MSstats(raw2)
  expect_true(all.equal(data.MSstats[, "precursorcharge"],
                        gsub(".*_([[:digit:]])$", "\\1", data.MSstats[, "fragmention"])))

  # test if negative values are replaced by NA
  raw2 <- raw
  raw2[raw2$fragmention == "1069078_FIIDPAAVITGR_2" &
       raw2$run == 3, "intensity"] <- -100
  data.MSstats <- convert_MSstats(raw2)
  expect_warning(convert_MSstats(raw2),
                 "Negative intensity values were replaced by NA")
  expect_true(is.na(data.MSstats[data.MSstats$fragmention == "1069078_FIIDPAAVITGR_2" &
                                 data.MSstats$run == 3, "intensity"]))
})

test_that("non-proteotypic peptides", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data$ProteinName <- gsub("1/Protein6", "2/Protein6/DECOY_ProteinX", data$ProteinName)
  data$ProteinName <- gsub("1/Protein5", "2/DECOY_ProteinX/Protein6", data$ProteinName)

  expect_true(sum(data$ProteinName == "2/Protein6/DECOY_ProteinX")>1)
  expect_true(sum(data$ProteinName == "2/DECOY_ProteinX/Protein6")>1)

  data <- unifyProteinGroupLabels(data)

  expect_false(sum(data$proteinname == "2/Protein6/DECOY_ProteinX")>1)
  expect_true(sum(data$proteinname == "2/DECOY_ProteinX/Protein6")>1)
  expect_false(sum(data$proteinname == "1/Protein6")>1)

  data <- removeDecoyProteins(data)

  expect_false(sum(data$proteinname == "2/Protein6/DECOY_ProteinX")>1)
  expect_false(sum(data$proteinname == "2/DECOY_ProteinX/Protein6")>1)
  expect_true(sum(data$proteinname == "1/Protein6")>1)
})
