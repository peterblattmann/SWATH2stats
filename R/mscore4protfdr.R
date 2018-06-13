#' Calculate the mscore which provides the given protein fdr hit rate.
#'
#' Hmmm the 3 functions in this family are largely equivalent, differing only in
#' the column which is used to query the name.  At some point in the near future
#' I would like to replace these with a generic function that can handle either
#' assays, peptides, or proteins.
#'
#' @param data SWATH2stats data to poke.
#' @param FFT I dunno.
#' @param fdr_target  What fdr to query?
#' @return an mscore which gets to this fdr hit rate.
#' @export
mscore4protfdr <- function(data, FFT=1, fdr_target=0.02) {
  # generate high resolution mscore levels to assess mscore cutoff for a given fdr_target
  mscore_levels_highres <- 10 ^ -(c(seq(2, 20, 0.05)))
  target.protein.highres <- NULL
  decoy.protein.highres <- NULL
  for (i in 1:length(mscore_levels_highres)) {
    target.protein.highres[i] <- length(
      unique(data[data[["decoy"]] == FALSE &
                  data[["m_score"]] <= mscore_levels_highres[i], c("proteinname")]))
    decoy.protein.highres[i] <- length(
      unique(data[data[["decoy"]] == TRUE &
                  data[["m_score"]] <= mscore_levels_highres[i], c("proteinname")]))
  }
  protein.fdr.highres <- (decoy.protein.highres / target.protein.highres) * FFT

  # pick mscore cutoff closest to (<=) fdr_target % protein FDR & filter data accordingly
  mscore_chosen <- mscore_levels_highres[protein.fdr.highres <= fdr_target][1]
  protein_fdr_chosen <- protein.fdr.highres[protein.fdr.highres <= fdr_target][1]
  message("Target protein FDR: ", fdr_target)
  message("Required overall m-score cutoff: ",
          signif(mscore_chosen, digits=5),"\n", "achieving protein FDR: ",
          signif(protein_fdr_chosen, digits=3))
  return(mscore_chosen)
}
