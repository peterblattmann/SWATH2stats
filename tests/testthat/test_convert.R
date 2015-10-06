context("data conversion")

test_that("data conversion", {
  data(OpenSWATH_data, package="SWATH2stats")
  data <- OpenSWATH_data
  data(Study_design, package="SWATH2stats")
  
  data <- sample_annotation(data, Study_design)
  data.filtered.mscore <- filter_mscore_requant(data, 0.01, 0.8)
  

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

  raw <- disaggregate(data.max)
  expect_true(subset(raw, FragmentIon == "1069078_FIIDPAAVITGR_2" & Run == 6)[,"Intensity"] == 188569)
  expect_true(subset(raw, FragmentIon == "279546_FPSIVGVAR_2" & Run == 5)[,"Intensity"] == 10805)
  
  
  data.MSstats <- convert4MSstats(raw)
  expect_that(dim(data.MSstats), equals(c(1800,10)))

  data.mapDIA <- convert4mapDIA(raw, RT = TRUE)
  expect_that(dim(data.mapDIA), equals(c(300,10)))
  
  data.aLFQ <- convert4aLFQ(raw)
  expect_that(dim(data.aLFQ), equals(c(1800,8)))
  
})

