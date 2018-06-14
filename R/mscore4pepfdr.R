#' Calculate the mscore which provides the given peptide fdr hit rate.
#'
#' @param data SWATH2stats data to poke.
#' @param FFT I dunno.
#' @param fdr_target  What fdr to query?
#' @return an mscore which gets to this fdr hit rate.
#' @export
mscore4pepfdr <- function(data, FFT=1, fdr_target=0.01) {
  ## generate high resolution mscore levels to assess mscore cutoff for a given fdr_target
  mscore_levels_highres <- 10 ^ -(c(seq(2, 20, 0.05)))
  target_peptides_highres <- NULL
  decoy_peptides_highres <- NULL
  for (i in 1:length(mscore_levels_highres)) {
    target_peptides.highres[i] <- length(
      unique(data[data[["decoy"]] == FALSE &
                  data[["m_score"]] <= mscore_levels_highres[i], c("fullpeptidename")]))
    decoy_peptides.highres[i] <- length(
      unique(data[data[["decoy"]] == TRUE &
                  data[["m_score"]] <= mscore_levels_highres[i], c("fullpeptidename")]))
  }
  peptide_fdr_highres <- (decoy_peptides_highres / target_peptides_highres) * FFT

  # pick mscore cutoff closest to (<=) fdr_target % peptide FDR & filter data accordingly
  mscore_chosen <- mscore_levels_highres[peptide_fdr_highres <= fdr_target][1]
  peptide_fdr_chosen <- peptide_fdr_highres[peptide_fdr_highres <= fdr_target][1]
  message("Target peptide FDR: ", fdr_target)
  message("Required overall m-score cutoff: ", signif(mscore_chosen, digits=5))
  message("Achieving peptide FDR: ", signif(peptide_fdr_chosen, digits=3))
  return(mscore_chosen)
}
