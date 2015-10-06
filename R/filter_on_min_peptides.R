filter_on_min_peptides <- function(data, n_peptides){
  data.prot.pep <- unique(data[,c("ProteinName", "FullPeptideName")])
  
  data.prot.pep.n <- tapply(data.prot.pep$FullPeptideName, data.prot.pep$ProteinName, length)
  
  
  prot.pep.names <- names(data.prot.pep.n[data.prot.pep.n >= n_peptides & !is.na(data.prot.pep.n)])
  data.2 <- data[data$ProteinName %in% prot.pep.names,]
  return(data.2)
}