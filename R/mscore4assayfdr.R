mscore4assayfdr <- function(data, FFT = 1, fdr_target = 0.01, mscore.col = "m_score") {
    mscore.col <- JPP_update(data, mscore.col)
    
    # generate high resolution mscore levels to assess mscore cutoff for a given
    # fdr_target
    mscore_levels_highres = 10^-(c(seq(2, 20, 0.05)))
    target.assays.highres <- NULL
    decoy.assays.highres <- NULL
    for (i in seq_len(length(mscore_levels_highres))) {
        target.assays.highres[i] <- length(unique(data[data$decoy == FALSE & data[, 
            mscore.col] <= mscore_levels_highres[i], c("transition_group_id")]))
        decoy.assays.highres[i] <- length(unique(data[data$decoy == TRUE & data[, 
            mscore.col] <= mscore_levels_highres[i], c("transition_group_id")]))
    }
    assay.fdr.highres <- (decoy.assays.highres/target.assays.highres) * FFT
    
    # pick mscore cutoff closest to (<=) fdr_target % peptide FDR & report
    mscore_chosen <- mscore_levels_highres[assay.fdr.highres <= fdr_target][1]
    assay_fdr_chosen <- assay.fdr.highres[assay.fdr.highres <= fdr_target][1]
    message("Target assay FDR: ", fdr_target, "\n")
    message("Required overall m-score cutoff:", signif(mscore_chosen, digits = 5), 
        "\n", "achieving assay FDR =", signif(assay_fdr_chosen, digits = 3), "\n")
    return(mscore_chosen)
}
