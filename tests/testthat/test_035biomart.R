context("biomart")

test_that("convert_ids", {
  data_table <- data.frame(Protein = c("Q01581", "P49327", "2/P63261/P60709"), 
                           Abundance = c(100, 3390, 43423))
  result <- convert_protein_ids(data_table)
  expect_that(result[which(result$Protein == "Q01581"), "hgnc_symbol"], equals("HMGCS1"))
  expect_that(result[which(result$Protein == "P49327"), "hgnc_symbol"], equals("FASN"))
 })
