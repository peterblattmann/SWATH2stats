utils::globalVariables(c("write.csv"))

assess_fdr_byrun <- function(data, FFT = 1, n.range = 20, output = "pdf_csv", plot = TRUE, filename = "FDR_report_byrun")
{
  # create m_score intervals to be tested
  test_levels <- 10^-seq(1:n.range)

  ## Identify the minimal m-score cutoff at which all runs still contain decoys
  decoy_count_lengths <- NULL
  for (i in seq(1:n.range)) {
    decoy_count_lengths[i] <- length(by(data[data$decoy == TRUE & data$m_score <= test_levels[i],
        c("transition_group_id")], data[data$decoy == TRUE & data$m_score <= test_levels[i], c("run_id")], length))
  }

  mscore_limit <- length(decoy_count_lengths[decoy_count_lengths == length(unique(data$run_id))])

  mscore_levels <- 10^-c(seq(2, mscore_limit))

  fdr_cube <- array(NA, dim = c(12, length(unique(data$run_id)), length(mscore_levels)))

  dimnames(fdr_cube) <- list(c("target_assays", "decoy_assays", "false_target_assays",
        "assay_fdr", "target_peptides", "decoy_peptides", "false_target_peptides", "peptide_fdr",
        "target_proteins", "decoy_proteins", "false_target_proteins", "protein_fdr"),
        sort(unique(data$run_id)), mscore_levels)

  length_unique <- function(X) {
    length(unique(X))
  }

  data.t <- data[data$decoy == 0, ]
  data.d <- data[data$decoy == 1, ]

  # for each m_score cutoff in mscore_levels, count targets and decoys and calculate false targets and fdr
  for (i in 1:length(mscore_levels)) {
    # for each run_id, count target (and decoy) assays identified ("id" column entries)
    # and store in pane i /row 1 (targets) /row 2 (decoys) & calculate false targets & fdr based on FFT
    fdr_cube[1, , i] <- by(data.t[data.t$m_score <= mscore_levels[i], c("transition_group_id")],
                           data.t[data.t$m_score <= mscore_levels[i],c("run_id")], length)
    fdr_cube[2, , i] <- by(data.d[data.d$m_score <= mscore_levels[i], c("transition_group_id")],
                           data.d[data.d$m_score <= mscore_levels[i], c("run_id")], length)
    fdr_cube[3, , i] <- fdr_cube[2, , i] * FFT
    fdr_cube[4, , i] <- fdr_cube[3, , i]/fdr_cube[1, , i]

    # for each run_id, count target (and decoy) peptides identified (unique "FullPeptideName" column entries)
    # and store in pane i /row 1 (targets) /row 2 (decoys) & calculate false targets & fdr based on FFT
    fdr_cube[5, , i] <- by(data.t[data.t$m_score <= mscore_levels[i], c("FullPeptideName")],
                           data.t[data.t$m_score <= mscore_levels[i], c("run_id")], length_unique)
    fdr_cube[6, , i] <- by(data.d[data.d$m_score <= mscore_levels[i], c("FullPeptideName")],
                           data.d[data.d$m_score <= mscore_levels[i], c("run_id")], length_unique)
    fdr_cube[7, , i] <- fdr_cube[6, , i] * FFT
    fdr_cube[8, , i] <- fdr_cube[7, , i]/fdr_cube[5, , i]

    # for each run_id, count target (and decoy) proteins identified (unique "ProteinName" column entries)
    # and store in pane i /row 1 (targets) /row 2 (decoys) & calculate false targets & fdr based on FFT
    fdr_cube[9, , i] <- by(data.t[data.t$m_score <= mscore_levels[i], c("ProteinName")],
                           data.t[data.t$m_score <= mscore_levels[i], c("run_id")], length_unique)
    fdr_cube[10, , i] <- by(data.d[data.d$m_score <= mscore_levels[i], c("ProteinName")],
                            data.d[data.d$m_score <= mscore_levels[i], c("run_id")], length_unique)
    fdr_cube[11, , i] <- fdr_cube[10, , i] * FFT
    fdr_cube[12, , i] <- fdr_cube[11, , i]/fdr_cube[9, , i]
  }

  # print fdr values to console (m_score cutoff 1e-2)
  message("The average FDR by run on assay level is ", round(mean(fdr_cube[4, , 1], na.rm=TRUE), digits=3), "\n")
  message("The average FDR by run on peptide level is ", round(mean(fdr_cube[8, , 1], na.rm=TRUE), digits=3), "\n")
  message("The average FDR by run on protein level is ", round(mean(fdr_cube[12, , 1], na.rm=TRUE), digits=3), "\n")

  if(output == "pdf_csv"){
    message("Individual run FDR qualities can be retrieved from ", paste(filename, ".csv"), "\n", sep="")
    # Write csv reports for mscore 1e-2 and 1e-3
    write.csv(fdr_cube[, , 1], file = paste(filename, "table_mscore_1e-2.csv", sep = "_"))
    write.csv(fdr_cube[, , 2], file = paste(filename, "table_mscore_1e-3.csv", sep = "_"))
    message(filename, ".csv reports written to working folder", "\n")
  }


  fdr_cube2 <- fdr_cube
  class(fdr_cube2) <- "fdr_cube"

  if(isTRUE(plot)){
    plot.fdr_cube(fdr_cube2, output = output, filename = filename)
  }

  if(output == "Rconsole"){
    return(fdr_cube2)
  }

}