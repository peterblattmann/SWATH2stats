mscore4pepfdr<- function(data, FFT = 1, fdr_target = 0.01, mscore.col = "m_score"){
  
  mscore.col <- JPP_update(data, mscore.col)
  
  # generate high resolution mscore levels to assess mscore cutoff for a given fdr_target
  mscore_levels_highres=10^-(c(seq(2, 20, 0.05)))
  target.peptides.highres<-NULL
  decoy.peptides.highres<-NULL
  for (i in seq_len(length(mscore_levels_highres)))
  {  
    target.peptides.highres[i]<-length(unique(data[data$decoy == FALSE & data[,mscore.col] <= mscore_levels_highres[i], c("FullPeptideName")]))
    decoy.peptides.highres[i]<-length(unique(data[data$decoy == TRUE & data[,mscore.col] <= mscore_levels_highres[i], c("FullPeptideName")]))
  }
  peptide.fdr.highres<-(decoy.peptides.highres/target.peptides.highres)*FFT
  
  # pick mscore cutoff closest to (<=) fdr_target % peptide FDR & filter data accordingly
  mscore_chosen<-mscore_levels_highres[peptide.fdr.highres<=fdr_target][1]
  peptide_fdr_chosen<-peptide.fdr.highres[peptide.fdr.highres<=fdr_target][1]
  message("Target peptide FDR:", fdr_target, "\n")
  message("Required overall m-score cutoff:", signif(mscore_chosen, digits=5), "\n", "achieving peptide FDR =", signif(peptide_fdr_chosen, digits=3), "\n")
  return(mscore_chosen)
}
