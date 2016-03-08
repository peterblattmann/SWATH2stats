utils::globalVariables(c("write.csv"))

assess_fdr_overall <-function(data, FFT = 1, n.range = 20, output="pdf_csv", plot=TRUE, filename="FDR_report_overall")
{
  mscore_levels<-10^-seq(n.range)

  # create vectors to store count results
  target.assays<-NULL
  decoy.assays<-NULL
  target.peptides<-NULL
  decoy.peptides<-NULL
  target.proteins<-NULL
  decoy.proteins<-NULL

  # loop over IDs with different m_score thresholds counting targets & decoys
  for (i in seq_len(length(mscore_levels)))
  {
    target.assays[i]<-length(unique(data[data$decoy == FALSE & data$m_score <= mscore_levels[i], c("transition_group_id")]))
    decoy.assays[i]<-length(unique(data[data$decoy == TRUE & data$m_score <= mscore_levels[i], c("transition_group_id")]))
    target.peptides[i]<-length(unique(data[data$decoy == FALSE & data$m_score <= mscore_levels[i], c("FullPeptideName")]))
    decoy.peptides[i]<-length(unique(data[data$decoy == TRUE & data$m_score <= mscore_levels[i], c("FullPeptideName")]))
    target.proteins[i]<-length(unique(data[data$decoy == FALSE & data$m_score <= mscore_levels[i], c("ProteinName")]))
    decoy.proteins[i]<-length(unique(data[data$decoy == TRUE & data$m_score <= mscore_levels[i], c("ProteinName")]))
  }

  # calculate false target fraction at cutoff (FDR) by decoy fraction * FFT
  assay.fdr<-(decoy.assays/target.assays)*FFT
  peptide.fdr<-(decoy.peptides/target.peptides)*FFT
  protein.fdr<-(decoy.proteins/target.proteins)*FFT
  # calculate estimated number of true identifications among the target hits
  true.target.assays<-target.assays-(decoy.assays*FFT)
  true.target.peptides<-target.peptides-(decoy.peptides*FFT)
  true.target.proteins<-target.proteins-(decoy.proteins*FFT)


  fdr_table <- structure(list(
    mscore_cutoff = mscore_levels,
    target.assays = target.assays,
    decoy.assays = decoy.assays,
    assay.fdr = assay.fdr,
    true.target.assays = true.target.assays,
    target.peptides = target.peptides,
    decoy.peptides = decoy.peptides,
    peptide.fdr = peptide.fdr,
    true.target.peptides = true.target.peptides,
    target.proteins = target.proteins,
    decoy.proteins = decoy.proteins,
    protein.fdr = protein.fdr,
    true.target.proteins = true.target.proteins),
    class="fdr_table")

  if(isTRUE(plot)){
    plot.fdr_table(fdr_table, filename= filename, output = output)
  }

  if(output == "Rconsole"){
    return(fdr_table)
  }

  if(output == "pdf_csv"){
    fdr_table_csv <- cbind(mscore_levels, target.assays, decoy.assays, assay.fdr, true.target.assays, target.peptides, decoy.peptides, peptide.fdr, true.target.peptides, target.proteins, decoy.proteins, protein.fdr, true.target.proteins)
    write.csv(fdr_table_csv, file=paste(filename, "table.csv", sep="_"), row.names=TRUE)
    message(filename,"table.csv written to working folder", "\n")
  }

}