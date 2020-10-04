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
#' @param upper_overall_peptide_fdr_limit  Option to relax or tighten the false
#'   discovery rate limit.
#' @param rm_decoy Logical T/F, whether decoy entries should be removed after
#'   the analysis. Defaults to TRUE. Can be useful to disable to track the
#'   influence on decoy fraction by further filtering steps such as requiring 2
#'   peptides per protein.
#' @param score_col Defines the column from which to retrieve the m_score. 
#'    If you use JPP (Rosenberger, Bludau et al. 2017) this can be used to 
#'    select between Protein and transition_group m_score.
#' @return Returns a data frame with the filtered data.
#' @author Moritz Heusel
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  data.fdr.filtered<-filter_mscore_fdr(data, FFT=0.7, 
#'                                       overall_protein_fdr_target=0.02,
#'                                       upper_overall_peptide_fdr_limit=0.1)
#' @export
filter_mscore_fdr <- function(data, 
                              FFT = 1, 
                              overall_protein_fdr_target = 0.02, 
                              upper_overall_peptide_fdr_limit = 0.05,
                              rm_decoy = TRUE, 
                              score_col = "m_score") {
    score_col <- JPP_update(data, score_col)

    mscore4protfdr_target <- mscore4protfdr(data, FFT, fdr_target = overall_protein_fdr_target)

    if (is.na(mscore4protfdr_target)) {
        stop("The overall_protein_fdr_target cannot be reached in this dataset.\n
         Consider using a higher accepted FDR criterion. \n
         Check mProphet models for target-decoy separation.")
    }

    # Initiate reporting to console
    message("filter_mscore_fdr is filtering the data...", "\n")
    message("-----------------------------------------------------------", "\n")
    message("finding m-score cutoff to achieve desired protein FDR in protein master list..",
        "\n")
    # Create master list at strict protein level FDR criterion
    protein_master_list <- unique(subset(data, data[, score_col] <= mscore4protfdr_target)$ProteinName)

    # Pre-Filter data based on upper_overall_peptide_fdr_limit
    message("finding m-score cutoff to achieve desired global peptide FDR..", "\n")
    data.f1 <- subset(data, data[, score_col] <= mscore4pepfdr(data, FFT, fdr_target = upper_overall_peptide_fdr_limit))

    # Filter prefiltered data down to entries mapping to the protein_master_list
    data.f2 <- subset(data.f1, data.f1$ProteinName %in% protein_master_list)

    # count remaining entries
    proteins <- length(protein_master_list)
    proteins.t <- length(unique(data.f2[data.f2$decoy == FALSE, c("ProteinName")]))
    proteins.d <- length(unique(data.f2[data.f2$decoy == TRUE, c("ProteinName")]))
    total.peptides <- length(unique(data.f1$FullPeptideName))
    total.peptides.t <- length(unique(data.f1[data.f2$decoy == FALSE, c("FullPeptideName")]))
    total.peptides.d <- length(unique(data.f1[data.f2$decoy == TRUE, c("FullPeptideName")]))
    mapping.peptides <- length(unique(data.f2$FullPeptideName))
    mapping.peptides.t <- length(unique(data.f2[data.f2$decoy == FALSE, c("FullPeptideName")]))
    mapping.peptides.d <- length(unique(data.f2[data.f2$decoy == TRUE, c("FullPeptideName")]))


    # print some numbers about the filtering results
    message("-------------------------------------------------------------", "\n")
    message("Proteins selected: ", "\n", "Total proteins selected: ", proteins, "\n",
        "Thereof target proteins: ", proteins.t, "\n", "Thereof decoy proteins: ",
        proteins.d, "\n")
    message("Peptides mapping to these protein entries selected:", "\n", "Total mapping peptides: ",
        mapping.peptides, "\n", "Thereof target peptides: ", mapping.peptides.t,
        "\n", "Thereof decoy peptides: ", mapping.peptides.d, "\n")
    message("Total peptides selected from:", "\n", "Total peptides: ", total.peptides,
        "\n", "Thereof target peptides: ", total.peptides.t, "\n", "Thereof decoy peptides: ",
        total.peptides.d, "\n")
    message("-------------------------------------------------------------", "\n")

    # test if all runs contain a decoy after peptide FDR filering in order to
    # calculate the local FDR
    n.run_id <- length(unique(data.f1$run_id))
    n.run_id_decoy <- length(unique(data.f1[data.f1$decoy == TRUE & data.f1[, score_col] <=
        0.01, ]$run_id))
    if (n.run_id == n.run_id_decoy) {

        fdr_cube <- assess_fdr_byrun(data.f1, FFT, output = "Rconsole", plot = FALSE)

        message("-------------------------------------------------------------",
            "\n")
        message("Individual run FDR quality of the peptides selected from:", "\n")
        message(signif(mean(fdr_cube[8, , 1]), 3))
    }
    if (n.run_id != n.run_id_decoy) {
        message("Individual run FDR quality of the peptides was not calculated",
            "\n", "as not every run contains a decoy.", "\n")
    }

    # return filtered data with or without decoy entries
    if (rm_decoy == FALSE) {
        message("The decoys have NOT been removed from the returned data", "\n")
        return(data.f2)
    } else {
        message("The decoys have been removed from the returned data", "\n")
        return(data.f2[data.f2$decoy == FALSE, ])
    }
}
