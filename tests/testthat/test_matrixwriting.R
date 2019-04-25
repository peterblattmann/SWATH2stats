context("Overview matrix writing")

test_that("Overview matrix writing", {
  data(OpenSWATH_data, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")
  data<-sample_annotation(OpenSWATH_data, Study_design)
  
  # peptide_matrix_test
  peptide_matrix <- write_matrix_peptides(data)
  
  expect_that(peptide_matrix[which(peptide_matrix$ProteinName_FullPeptideName == "1/iRT_protein_ATFGVDESNAEVK"),
                             "Hela_Control_1_1"],
              equals(86041))

  # test rm.decoy option
  peptide_matrix2 <- write_matrix_peptides(data, rm.decoy=TRUE)
  expect_true(length(grep("DECOY", peptide_matrix$ProteinName_FullPeptideName)) > 0)
  expect_false(length(grep("DECOY", peptide_matrix2$ProteinName_FullPeptideName)) > 0)
  
  # test write.csv option
  expect_message(write_matrix_peptides(data, write.csv = TRUE),
                 "Peptide overview matrix SWATH2stats_overview_matrix_peptidelevel.csv written to working folder.")
  
  
  # protein_matrix_test
  protein_matrix <- write_matrix_proteins(data)

  expect_that(protein_matrix[which(protein_matrix$ProteinName == "1/iRT_protein"),
                                "Hela_Control_1_1"], equals(1538879))

  # test rm.decoy option
  protein_matrix2 <- write_matrix_proteins(data, rm.decoy=TRUE)
  expect_true(length(grep("DECOY", protein_matrix$ProteinName)) > 0)
  expect_false(length(grep("DECOY", protein_matrix2$ProteinName)) > 0)
  
  # test write.csv option
  expect_message(write_matrix_proteins(data, write.csv = TRUE),
                 "Protein overview matrix SWATH2stats_overview_matrix_proteinlevel.csv written to working folder")
})
