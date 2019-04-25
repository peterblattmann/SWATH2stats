utils::globalVariables(c("aggregate"))

count_analytes <- function(data, column.levels = c("transition_group_id", "FullPeptideName", 
    "ProteinName"), column.by = "run_id", rm.decoy = TRUE) {
    if (sum(colnames(data) == "decoy") == 1 & isTRUE(rm.decoy)) {
        # data <- subset(data, decoy == 0)
        data <- data[data$decoy == 0, ]
    }
    
    data.n <- aggregate(data[, column.levels], by = list(data[, column.by]), function(x) length(unique(x)))
    colnames(data.n)[1] <- column.by
    return(data.n)
}
