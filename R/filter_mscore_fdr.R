filter_mscore_fdr <- function(data, FFT=1, overall_protein_fdr_target=0.02, mscore_limit=0.01,
                              upper_overall_peptide_fdr_limit=0.05, rm.decoy=TRUE) {
  mscore4protfdr_target <- mscore4protfdr(data, FFT, fdr_target=overall_protein_fdr_target)

  if (is.na(mscore4protfdr_target)) {
    stop("The overall_protein_fdr_target cannot be reached in this dataset. \n
         Consider using a higher accepted FDR criterion. \n
         Check mProphet models for target-decoy separation.")
  }

  # Initiate reporting to console
  message("filter_mscore_fdr is filtering the data...")
  message("finding m-score cutoff to achieve desired protein FDR in protein master list..")
  ## Create master list at strict protein level FDR criterion
  lower_idx <- data[["m_score"]] <= mscore4protfdr_target
  if (sum(lower_idx) == 0) {
    warning("No m scores were lower than the fdr target, that cannot be good.")
  }
  lower_fdr <- data[lower_idx, ]
  protein_master_list <- unique(lower_fdr[["proteinname"]])

  # Pre-Filter data based on upper_overall_peptide_fdr_limit
  message("finding m-score cutoff to achieve desired global peptide FDR..")
  lower_limit <- mscore4pepfdr(data, FFT, fdr_target=upper_overall_peptide_fdr_limit)
  lower_idx <- data[["m_score"]] <= lower_limit
  data.f1 <- data[lower_idx, ]

  ## Filter prefiltered data down to entries mapping to the protein_master_list
  found_idx <- data.f1[["proteinname"]] %in% protein_master_list
  data.f2 <- data.f1[found_idx, ]

  # count remaining entries
  proteins <- length(protein_master_list)
  proteins.t <- length(unique(data.f2[data.f2[["decoy"]] == FALSE, c("proteinname")]))
  proteins.d <- length(unique(data.f2[data.f2[["decoy"]] == TRUE, c("proteinname")]))
  total.peptides <- length(unique(data.f1[["fullpeptidename"]]))
  total.peptides.t <- length(unique(data.f1[data.f2[["decoy"]] == FALSE, c("fullpeptidename")]))
  total.peptides.d <- length(unique(data.f1[data.f2[["decoy"]] == TRUE, c("fullpeptidename")]))
  mapping.peptides <- length(unique(data.f2[["fullpeptidename"]]))
  mapping.peptides.t <- length(unique(data.f2[data.f2[["decoy"]] == FALSE, c("fullpeptidename")]))
  mapping.peptides.d <- length(unique(data.f2[data.f2[["decoy"]] == TRUE, c("fullpeptidename")]))


  # print some numbers about the filtering results
  message("Proteins selected: ", "\n",
      "Total proteins selected: ", proteins, "\n",
      "Final target proteins: ", proteins.t, "\n",
      "Final decoy proteins: ", proteins.d, "\n")
  message("Peptides mapping to these protein entries selected:", "\n",
      "Total mapping peptides: ", mapping.peptides, "\n",
      "Final target peptides: ", mapping.peptides.t, "\n",
      "Final decoy peptides: ", mapping.peptides.d, "\n")
  message("Total peptides selected from:", "\n",
      "Total peptides: ", total.peptides, "\n",
      "Final target peptides: ", total.peptides.t, "\n",
      "Final decoy peptides: ", total.peptides.d, "\n")

  # test if all runs contain a decoy after peptide FDR filering in order to calculate the local FDR
  n.run_id <- length(unique(data.f1[["run_id"]]))
  n.run_id_decoy <- length(
    unique(data.f1[data.f1[["decoy"]] == TRUE &
                   data.f1[["m_score"]] <= mscore_limit, ][["run_id"]]))
  if (n.run_id == n.run_id_decoy) {
    fdr_cube <- assess_fdr_byrun(data.f1, FFT, output="Rconsole", plot=FALSE)
    message("Individual run FDR quality of the peptides selected from:")
    message(signif(mean(fdr_cube[8, , 1]), 3))
  }
  if (n.run_id != n.run_id_decoy) {
    message("Individual run FDR quality of the peptides was not calculated",
            "\n", "as not every run contains a decoy.")
  }

  # return filtered data with or without decoy entries
  if (rm.decoy == FALSE) {
    message("The decoys have NOT been removed from the returned data.")
    return(data.f2)
  } else {
    message("The decoys have been removed from the returned data.")
    return(data.f2[data.f2[["decoy"]] == FALSE, ])
  }
}
