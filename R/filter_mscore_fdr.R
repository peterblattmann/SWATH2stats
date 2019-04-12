#' Filter annotated OpenSWATH/pyProphet output table to achieve a high FDR
#' quality data matrix with controlled overall protein FDR and quantitative
#' values for all peptides mapping to these high-confidence proteins (up to a
#' desired overall peptide level FDR quality).
#'
#' This function controls the protein FDR over a multi-run OpenSWATH/pyProphet
#' output table and filters all quantitative values to a desired overall/global
#' peptide FDR level.
#' It first finds a suitable m-score cutoff to minimally achieve a desired
#' global FDR quality on a protein master list based on the function
#' mscore4protfdr. It then finds a suitable m-score cutoff to minimally achieve
#' a desired global FDR quality on peptide level based on the function
#' mscore4pepfdr. Finally, it reports all the peptide quantities derived based
#' on the peptide level cutoff for only those peptides mapping to the protein
#' master list. It further summarizes the protein and peptide numbers remaining
#' after the filtering and evaluates the individual run FDR qualities of the
#' peptides (and quantitation events) selected.
#'
#' @param data  Annotated OpenSWATH/pyProphet data table.
#' @param FFT  Ratio of false positives to true negatives, q-values from
#'   [Injection_name]_full_stat.csv in pyProphet stats output. As an
#'   approximation, the q-values of multiple runs are averaged and supplied as
#'   argument FFT. Numeric from 0 to 1. Defaults to 1, the most conservative
#'   value (1 Decoy indicates 1 False target). For further details see the
#'   Vignette Section 1.3 and 4.1.
#' @param overall_protein_fdr_target FDR target for the protein master list for
#'   which quantitative values down to the less strict peptide_fdr criterion
#'   will be kept/reported. Defaults to 0.02.
#' @param mscore_limit  FDR target for the quantitative values kept/reported for
#'   all peptides mapping to the high-confidence protein master list. Defaults
#'   to 0.05. If all values up to m_score 0.01 shall be kept, set = 1.
#' @param upper_overall_peptide_fdr_limit  Option to relax or tighten the false
#'   discovery rate limit.
#' @param rm.decoy Logical T/F, whether decoy entries should be removed after
#'   the analysis. Defaults to TRUE. Can be useful to disable to track the
#'   influence on decoy fraction by further filtering steps such as requiring 2
#'   peptides per protein.
#' @return Returns a data frame with the filtered data.
#' @author Moritz Heusel
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.fdr.filtered<-filter_mscore_fdr(data, FFT=0.7, overall_protein_fdr_target=0.02,
#'                                       upper_overall_peptide_fdr_limit=0.1)
#' @export
filter_mscore_fdr <- function(data, FFT=1, overall_protein_fdr_target=0.02, mscore_limit=0.01,
                              upper_overall_peptide_fdr_limit=0.05, rm.decoy=TRUE,
                              mscore.col="m_score") {
  mscore.col <- JPP_update(data, mscore.col)
  mscore4protfdr_target <- mscore4protfdr(data, FFT, fdr_target=overall_protein_fdr_target)

  if (is.na(mscore4protfdr_target)) {
    stop("The overall_protein_fdr_target cannot be reached in this dataset. \n
         Consider using a higher accepted FDR criterion. \n
         Check mProphet models for target-decoy separation.")
  }

  ## Initiate reporting to console
  message("filter_mscore_fdr is filtering the data...")
  message("finding m-score cutoff to achieve desired protein FDR in protein master list..")
  ## Create master list at strict protein level FDR criterion
  lower_idx <- data[[mscore.col]] <= mscore4protfdr_target
  if (sum(lower_idx) == 0) {
    warning("No m scores were lower than the fdr target, that cannot be good.")
  }
  lower_fdr <- data[lower_idx, ]
  protein_master_list <- unique(lower_fdr[["proteinname"]])

  ## Pre-Filter data based on upper_overall_peptide_fdr_limit
  message("finding m-score cutoff to achieve desired global peptide FDR..")
  lower_limit <- mscore4pepfdr(data, FFT, fdr_target=upper_overall_peptide_fdr_limit)
  lower_idx <- data[[mscore.col]] <= lower_limit
  data.f1 <- data[lower_idx, ]

  ## Filter prefiltered data down to entries mapping to the protein_master_list
  found_idx <- data.f1[["proteinname"]] %in% protein_master_list
  data.f2 <- data.f1[found_idx, ]

  ## count remaining entries
  proteins <- length(protein_master_list)
  proteins.t <- length(unique(data.f2[data.f2[["decoy"]] == FALSE, c("proteinname")]))
  proteins.d <- length(unique(data.f2[data.f2[["decoy"]] == TRUE, c("proteinname")]))
  total.peptides <- length(unique(data.f1[["fullpeptidename"]]))
  total.peptides.t <- length(unique(data.f1[data.f2[["decoy"]] == FALSE, c("fullpeptidename")]))
  total.peptides.d <- length(unique(data.f1[data.f2[["decoy"]] == TRUE, c("fullpeptidename")]))
  mapping.peptides <- length(unique(data.f2[["fullpeptidename"]]))
  mapping.peptides.t <- length(unique(data.f2[data.f2[["decoy"]] == FALSE, c("fullpeptidename")]))
  mapping.peptides.d <- length(unique(data.f2[data.f2[["decoy"]] == TRUE, c("fullpeptidename")]))

  ## print some numbers about the filtering results
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

  ## test if all runs contain a decoy after peptide FDR filering in order to calculate the local FDR
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

  ## return filtered data with or without decoy entries
  if (rm.decoy == FALSE) {
    message("The decoys have NOT been removed from the returned data.")
    return(data.f2)
  } else {
    message("The decoys have been removed from the returned data.")
    return(data.f2[data.f2[["decoy"]] == FALSE, ])
  }
}
