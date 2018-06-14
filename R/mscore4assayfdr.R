#' Calculate the mscore which provides the given assay fdr hit rate.
#'
#' @param data SWATH2stats data to poke.
#' @param FFT I dunno.
#' @param fdr_target  What fdr to query?
#' @return an mscore which gets to this fdr hit rate.
#' @export
mscore4assayfdr <- function(data, FFT=1, fdr_target=0.01) {
  ## generate high resolution mscore levels to assess mscore cutoff for a given fdr_target
  mscore_levels_highres <- 10 ^ -(c(seq(2, 20, 0.05)))
  target_assays_highres <- NULL
  decoy_assays_highres <- NULL
  for (i in 1:length(mscore_levels_highres)) {
    target_assays_highres[i] <- length(unique(data[data[["decoy"]] == FALSE &
                                                   data[["m_score"]] <= mscore_levels_highres[i],
                                                   c("transition_group_id")]))
    decoy_assays_highres[i] <- length(unique(data[data[["decoy"]] == TRUE &
                                                  data[["m_score"]] <= mscore_levels_highres[i], c("transition_group_id")]))
  }
  assay_fdr_highres <- (decoy_assays_highres / target_assays_highres) * FFT

  # pick mscore cutoff closest to (<=) fdr_target % peptide FDR & report
  mscore_chosen <- mscore_levels_highres[assay_fdr_highres <= fdr_target][1]
  assay_fdr_chosen <- assay_fdr_highres[assay_fdr_highres <= fdr_target][1]
  message("Target assay FDR: ", fdr_target)
  message("Required overall m-score cutoff: ",
          signif(mscore_chosen, digits=5), "\n", "achieving assay FDR: ",
          signif(assay_fdr_chosen, digits=3))
  return(mscore_chosen)
}
