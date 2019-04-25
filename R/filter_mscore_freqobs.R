utils::globalVariables(c("Peptide_Charge", ".N"))

filter_mscore_freqobs <- function(data, mscore, percentage = NULL, rm.decoy = TRUE, 
    mscore.col = "m_score") {
    
    mscore.col <- JPP_update(data, mscore.col)
    
    data$Peptide_Charge <- paste(data$FullPeptideName, data$Charge)
    
    if (sum(colnames(data) == "decoy") == 1 & isTRUE(rm.decoy)) {
        # data <- subset(data, decoy == 0)
        data <- data[data$decoy == 0, ]
    }
    
    # data.filtered <- subset(data, m_score <= mscore)
    data.filtered <- data[data[, mscore.col] <= mscore, ]
    data.filtered <- data.table(data.filtered)
    
    data.filtered <- data.filtered[, c("Peptide_Charge", "aggr_Peak_Area"), with = FALSE]
    setkey(data.filtered, Peptide_Charge)
    data.n <- data.filtered[, .N, by = "Peptide_Charge"]
    
    if (is.null(percentage)) {
        percentage <- 0
    }
    
    threshold <- nlevels(factor(data$filename)) * percentage
    message("Treshold, peptides need to have been quantified in more conditions than: ", 
        threshold)
    
    peptides.filtered <- data.n[data.n$N >= threshold]
    peptides.filtered <- data.frame(Peptides_Charge = peptides.filtered$Peptide_Charge)
    
    message("Fraction of peptides selected: ", signif(length(unique(peptides.filtered$Peptides_Charge))/length(unique(data$Peptide_Charge)), 
        digits = 2))
    
    data.filtered <- merge(data, peptides.filtered, by.x = "Peptide_Charge", by.y = "Peptides_Charge")
    
    message("Dimension difference: ", paste(dim(data) - dim(data.filtered), collapse = ", "))
    
    return(data.filtered)
}

