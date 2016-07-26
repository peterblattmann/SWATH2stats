filter_mscore_fdr <- function(data, FFT = 1, overall_protein_fdr_target = 0.02, upper_overall_peptide_fdr_limit = 0.05, rm.decoy = TRUE)
{
  mscore4protfdr_target <- mscore4protfdr(data, FFT, fdr_target = overall_protein_fdr_target)
  
  if(is.na(mscore4protfdr_target)){
    stop("The overall_protein_fdr_target cannot be reached in this dataset. \n
         Consider using a higher accepted FDR criterion. \n
         Check mProphet models for target-decoy separation.")
  }
  
  # Initiate reporting to console
  message("filter_mscore_fdr is filtering the data...", "\n")
  message("-------------------------------------------------------------", "\n")
  message("finding m-score cutoff to achieve desired protein FDR in protein master list..", "\n")
    # Create master list at strict protein level FDR criterion
  protein_master_list <- unique(subset(data, data$m_score <= mscore4protfdr_target)$ProteinName)
  
  # Pre-Filter data based on upper_overall_peptide_fdr_limit
  message("finding m-score cutoff to achieve desired global peptide FDR..", "\n", "")
    data.f1 <- subset(data, data$m_score <=
    mscore4pepfdr(data, FFT, fdr_target = upper_overall_peptide_fdr_limit))

  # Filter prefiltered data down to entries mapping to the protein_master_list
  data.f2 <- subset(data.f1, data.f1$ProteinName %in% protein_master_list)

  # count remaining entries
  proteins<-length(protein_master_list)
  proteins.t<-length(unique(data.f2[data.f2$decoy == FALSE, c("ProteinName")]))
  proteins.d<-length(unique(data.f2[data.f2$decoy == TRUE, c("ProteinName")]))
  total.peptides<-length(unique(data.f1$FullPeptideName))
  total.peptides.t<-length(unique(data.f1[data.f2$decoy == FALSE, c("FullPeptideName")]))
  total.peptides.d<-length(unique(data.f1[data.f2$decoy == TRUE, c("FullPeptideName")]))
  mapping.peptides<-length(unique(data.f2$FullPeptideName))
  mapping.peptides.t<-length(unique(data.f2[data.f2$decoy == FALSE, c("FullPeptideName")]))
  mapping.peptides.d<-length(unique(data.f2[data.f2$decoy == TRUE, c("FullPeptideName")]))

  fdr_cube <- assess_fdr_byrun(data.f1, FFT, output = "Rconsole", plot = FALSE)

  # print some numbers about the filtering results
  message("-------------------------------------------------------------", "\n")
  message("Proteins selected: ", "\n",
      "Total proteins selected: ", proteins, "\n",
      "Thereof target proteins: ", proteins.t, "\n",
      "Thereof decoy proteins: ", proteins.d, "\n")
  message("Peptides mapping to these protein entries selected:", "\n",
      "Total mapping peptides: ", mapping.peptides, "\n",
      "Thereof target peptides: ", mapping.peptides.t, "\n",
      "Thereof decoy peptides: ", mapping.peptides.d, "\n")
  message("Total peptides selected from:", "\n",
      "Total peptides: ", total.peptides, "\n",
      "Thereof target peptides: ", total.peptides.t, "\n",
      "Thereof decoy peptides: ", total.peptides.d, "\n")
  message("-------------------------------------------------------------", "\n")
  
  # test if all runs contain a decoy after peptide FDR filering in order to calculate the local FDR
  n.run_id <- length(unique(data.f1$run_id)) 
  n.run_id_decoy <- length(unique(data.f1[data.f1$decoy == TRUE & data.f1$m_score <= 0.01,]$run_id))
  if(n.run_id == n.run_id_decoy){
    fdr_cube <- assess_fdr_byrun(data.f1, FFT, output = "Rconsole", 
                                 plot = FALSE)
    
    message("-------------------------------------------------------------", 
            "\n")
    message("Individual run FDR quality of the peptides selected from:", 
            "\n")
    message(signif(mean(fdr_cube[8, , 1]), 3))
  }
  if(n.run_id != n.run_id_decoy){
    message("Individual run FDR quality of the peptides was not calculated", 
            "\n", "as not every run contains a decoy.", 
            "\n")
  }

  # return filtered data with or without decoy entries
  if (rm.decoy == FALSE){
    message("The decoys have NOT been removed from the returned data", "\n")
    return(data.f2)
  }
  else{
    message("The decoys have been removed from the returned data", "\n")
    return(data.f2[data.f2$decoy == FALSE,])
  }
}