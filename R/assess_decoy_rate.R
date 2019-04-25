assess_decoy_rate <- function(data) {
    if (sum(colnames(data) == "decoy") < 1) {
        stop("There is no decoy column in the table")
    }
    if (sum(colnames(data) == "FullPeptideName") < 1) {
        stop("There is no FullPeptideName column in the table")
    }
    
    add.colnames <- colnames(data)
    add.colnames <- add.colnames[add.colnames != "decoy"]
    
    .non_decoy.peptides <- unique(data[data$decoy == FALSE, c("FullPeptideName")])
    .decoy.peptides <- unique(data[data$decoy == TRUE, c("FullPeptideName")])
    
    message("Number of non-decoy peptides: ", length(.non_decoy.peptides), "\n", 
        "Number of decoy peptides: ", length(.decoy.peptides), "\n", "Decoy rate: ", 
        sprintf("%.4f", (length(.decoy.peptides)/length(.non_decoy.peptides))))
}
