context("data conversion")

test_that("data conversion", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data(Study_design, package="SWATH2stats")

  data <- sample_annotation(data, Study_design)
  data.filtered.mscore <- filter_mscore_freqobs(data, 0.01, 0.8)

  # test if gives warning when data is in data.table format
  data2 <- data.table(data)
  expect_error(sample_annotation(data2, Study_design))

  data.proteotypic <- filter_proteotypic_peptides(data.filtered.mscore)
  expect_that(filter_proteotypic_peptides(data.filtered.mscore), shows_message("Number of proteins detected: 11"))
  expect_that(filter_proteotypic_peptides(data.filtered.mscore), shows_message("Number of proteotypic peptides detected: 131"))
  expect_that(dim(data.proteotypic), equals(c(972,62)))


  data.all <- filter_all_peptides(data.filtered.mscore)
  expect_that(filter_proteotypic_peptides(data.filtered.mscore), shows_message("Number of proteins detected: 11"))
  expect_that(dim(data.all), equals(c(984,62)))

  data.max <- filter_on_max_peptides(data.filtered.mscore, 5)
  expect_that(dim(data.max), equals(c(300,62)))
  expect_that(filter_on_max_peptides(data.filtered.mscore, 5), shows_message("Before filtering: \n  Number of proteins: 10\n  Number of peptides: 133\n\nPercentage of peptides removed: 69.17%\n\nAfter filtering: \n  Number of proteins: 10\n  Number of peptides: 41\n\n"))

  data.min <- filter_on_min_peptides(data.filtered.mscore, 5)
  expect_that(dim(data.min), equals(c(918,62)))

  data.python <- convert4pythonscript(data.max)
  expect_that(identical(grep("UniMod\\:", data.max[,"aggr_Fragment_Annotation"]), grep("UniMod\\_", data.python[,"aggr_Fragment_Annotation"])), is_true())
  expect_that(identical(grep("UniMod\\:", data.max[,"aggr_Fragment_Annotation"]), grep("UniMod\\_", data.max[,"aggr_Fragment_Annotation"])), is_false())
  expect_that(identical(grep("UniMod\\:", data.max[,"FullPeptideName"]), grep("UniMod\\_", data.python[,"FullPeptideName"])), is_true())
  expect_that(identical(grep("UniMod\\:", data.max[,"FullPeptideName"]), grep("UniMod\\_", data.max[,"FullPeptideName"])), is_false())

  raw <- disaggregate(data.max)
  expect_true(subset(raw, FragmentIon == "1069078_FIIDPAAVITGR_2" & Run == 6)[,"Intensity"] == 188569)
  expect_true(subset(raw, FragmentIon == "279546_FPSIVGVAR_2" & Run == 5)[,"Intensity"] == 10805)

  # test error if transitions and data don't match
  data.max.test1 <- data.max
  data.max.test1$aggr_Peak_Area <- as.character(data.max.test1$aggr_Peak_Area)
  data.max.test1[data.max.test1$Run == 6 & data.max.test1$BioReplicate == 3 & data.max.test1$Peptide_Charge == "FIIDPAAVITGR 2","aggr_Peak_Area"] <- "188569.000000;81947.000000;77295.000000;50218.000000"
  expect_error(disaggregate(data.max.test1))
  #expect_that(disaggregate(data.max.test1), throws_error("Error in disaggregate(data.max.test1) : \n  The number of transitions annotated and measured do not match in the following transitions:\n 1069078_FIIDPAAVITGR_2, 1069082_FIIDPAAVITGR_2, 1069080_FIIDPAAVITGR_2, 1069071_FIIDPAAVITGR_2, 1069084_FIIDPAAVITGR_2, 1069076_FIIDPAAVITGR_2\n", fixed=TRUE))
  expect_that(disaggregate(data.max.test1), throws_error("The number of transitions annotated and measured do not match in the following transitions:\n 1069078_FIIDPAAVITGR_2, 1069082_FIIDPAAVITGR_2, 1069080_FIIDPAAVITGR_2, 1069071_FIIDPAAVITGR_2, 1069084_FIIDPAAVITGR_2, 1069076_FIIDPAAVITGR_2", fixed=FALSE))
  
  # test if spectral library does not contain same number of transitions
  data.max.test2 <- data.max.test1
  data.max.test2$aggr_Fragment_Annotation <- as.character(data.max.test2$aggr_Fragment_Annotation)
  data.max.test2[data.max.test1$Run == 6 & data.max.test1$BioReplicate == 3 & data.max.test1$Peptide_Charge == "FIIDPAAVITGR 2","aggr_Fragment_Annotation"] <- "1069078_FIIDPAAVITGR_2;1069082_FIIDPAAVITGR_2;1069080_FIIDPAAVITGR_2;1069071_FIIDPAAVITGR_2"
  expect_that(disaggregate(data.max.test2), shows_message("The library contains between 4 and 6 transitions per precursor.", fixed=TRUE))
  expect_that(dim(disaggregate(data.max.test2)), equals(c(1798,10)))

  data.MSstats <- convert4MSstats(raw)
  expect_that(dim(data.MSstats), equals(c(1800,10)))

  data.mapDIA <- convert4mapDIA(raw, RT = TRUE)
  expect_that(dim(data.mapDIA), equals(c(300,10)))

  data.aLFQ <- convert4aLFQ(raw)
  expect_that(dim(data.aLFQ), equals(c(1800,8)))
  
  expect_warning(convert4aLFQ(raw[c(1:6),]), "The aLFQ package should only be used with transition-level data.")
  
  # test if warning is displayed when there are several values for a data point
  raw2 <- raw
  raw2$Intensity <- raw$Intensity + 20
  raw2 <- rbind(raw[raw$FragmentIon == "1069078_FIIDPAAVITGR_2" & raw$Condition %in% c("Hela_Treatment"),], raw2)
  expect_warning(convert4mapDIA(raw2, RT = TRUE), "Data contains several intensity values per condition")

  # test if warning is displayed when column is missing
  raw2 <- raw
  raw2$RT <- NULL
  expect_warning(convert4mapDIA(raw2, RT = TRUE), 
                 "One or several columns required by mDIA were not in the data and filled with NAs")
  
  # test if extraction of Precursor charge works
  raw2 <- raw
  raw2$PrecursorCharge <- NULL
  expect_message(convert4MSstats(raw2), 
                 "Missing columns: PrecursorCharge")
  data.MSstats <- convert4MSstats(raw2)
  expect_true(all.equal(data.MSstats[,"PrecursorCharge"], gsub(".*_([[:digit:]])$", "\\1", data.MSstats[, "FragmentIon"])))
  
  # test if negative values are replaced by NA
  raw2 <- raw
  raw2[raw2$FragmentIon == "1069078_FIIDPAAVITGR_2" & raw2$Run == 3, "Intensity"] <- -100
  data.MSstats <- convert4MSstats(raw2)
  expect_warning(convert4MSstats(raw2), 
                 "Negative intensity values were replaced by NA")
  expect_true(is.na(data.MSstats[data.MSstats$FragmentIon == "1069078_FIIDPAAVITGR_2" & data.MSstats$Run == 3, "Intensity"]))
  

})


test_that("non-proteotypic peptides", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data$ProteinName <- gsub("1/Protein6", "2/Protein6/DECOY_ProteinX", data$ProteinName)
  data$ProteinName <- gsub("1/Protein5", "2/DECOY_ProteinX/Protein6", data$ProteinName)
  
  expect_true(sum(data$ProteinName == "2/Protein6/DECOY_ProteinX")>1)
  expect_true(sum(data$ProteinName == "2/DECOY_ProteinX/Protein6")>1)
  
  data <- unifyProteinGroupLabels(data)
  
  expect_false(sum(data$ProteinName == "2/Protein6/DECOY_ProteinX")>1)
  expect_true(sum(data$ProteinName == "2/DECOY_ProteinX/Protein6")>1)
  expect_false(sum(data$ProteinName == "1/Protein6")>1)
  
  data <- removeDecoyProteins(data)
  
  expect_false(sum(data$ProteinName == "2/Protein6/DECOY_ProteinX")>1)
  expect_false(sum(data$ProteinName == "2/DECOY_ProteinX/Protein6")>1)
  expect_true(sum(data$ProteinName == "1/Protein6")>1)
})
