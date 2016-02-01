context("Overview matrix writing")

test_that("Overview matrix writing", {
  data(OpenSWATH_data, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")
  data<-sample_annotation(OpenSWATH_data, Study_design)
  peptide_matrix <- write_matrix_peptides(data)
  # peptide_matrix_test
  expect_that(peptide_matrix[which(peptide_matrix$ProteinName_FullPeptideName == "1/iRT_protein_ATFGVDESNAEVK"),
                             "Hela_Control_1_1"],
              equals(86041))

  # assess_fdr_overall output test
  protein_matrix <- write_matrix_proteins(data)

  expect_that(protein_matrix[which(protein_matrix$ProteinName == "1/iRT_protein"),
                                "Hela_Control_1_1"], equals(1538879))

})
