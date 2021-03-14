utils::globalVariables(c("Var1", "Var2", "value"))

#' Plots the correlation between injections.
#'
#' This function plots the Pearson's and Spearman correlation between
#' samples. If decoys are present these are removed before plotting.
#'
#' @param data Data frame that is produced by the OpenSWATH/pyProphet workflow.
#' @param column_values   Indicates the columns for which the correlation is
#'   assessed. This can be the Intensity or Signal, but also the retention
#'   time.
#' @param comparison  The comparison for assessing the variability. Default is
#'   to assess the variability per transition_group_id over the different
#'   Condition and Replicates. Comparison is performed using the dcast()
#'   function of the reshape2 package.
#' @param fun_aggregate  If for the comparison values have to be aggregated one
#'   needs to provide the function here.
#' @param label  Option to print correlation value in the plot.
#' @param ... Further arguments passed to methods.
#' @return Plots in Rconsole a correlation heatmap and returns the data frame
#'   used to do the plotting.
#' @author Peter Blattmann
#' @examples{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  information <- plot_correlation_between_samples(data)
#' }
#' @importFrom reshape2 dcast melt
#' @importFrom ggplot2 aes geom_tile scale_fill_gradient xlab ylab labs 
#'  geom_text theme theme_bw scale_x_discrete scale_y_discrete
#' @importFrom stats cor
#' @export
plot_correlation_between_samples <- function(data, 
                                             column_values = "Intensity",
                                             comparison = transition_group_id ~ Condition + BioReplicate,
                                             fun_aggregate = NULL, 
                                             label = TRUE, ...) {
    if (sum(colnames(data) == "decoy") == 1) {
        data <- data[data$decoy == 0, ]
    }
    data.c <- dcast(data, comparison, value.var = column_values, fun.aggregate = fun_aggregate)
    dep.vars <- length(all.vars(comparison[[2]]))
    indep.vars <- length(all.vars(comparison[[3]]))
    
    pearson.cor <- cor(data.c[, (dep.vars + 1):dim(data.c)[2]], 
                       use = "pairwise.complete.obs", method = "pearson")
    pearson.cor[lower.tri(pearson.cor)] <- NA
    
    spearman.cor <- cor(data.c[, (dep.vars + 1):dim(data.c)[2]], 
                        use = "pairwise.complete.obs", method = "spearman")
    spearman.cor[upper.tri(spearman.cor, diag = TRUE)] <- NA
    
    pearson.cor <- melt(pearson.cor)
    pearson.cor$method <- "pearson"
    spearman.cor <- melt(spearman.cor)
    spearman.cor$method <- "spearman"
    
    data.plot <- rbind(pearson.cor, spearman.cor)
    data.plot <- data.plot[!is.na(data.plot$value), ]
    
    if (label) {
        p <- (ggplot(data.plot, aes(x = Var2, y = Var1, fill = value)) +
                  geom_tile() +
                  scale_fill_gradient(low = "white", high = "red",
                                      name = "Correlation\n[R or rho]") +
                  xlab("") +
                  ylab("") +
                  labs(title = paste0(column_values, 
                                      " correlation between samples:\nPearson (upper triangle) and Spearman correlation (lower triangle)")) +
                  geom_text(aes_string(label = "round(data.plot$value, digits = 2)")) +
                  scale_x_discrete(expand = c(0, 0)) +
                  scale_y_discrete(limits = rev(levels(data.plot$Var1)),
                                   expand = c(0, 0)) +
                  theme_bw() +
                  theme(plot.title = element_text(hjust = 0.5, vjust = 1),
                        axis.text.x = element_text(angle = 90, hjust = 1)))
    } else {
        p <- (ggplot(data.plot, aes_string(x = "Var2", y = "Var1", fill = "value")) +
                  geom_tile() +
                  scale_fill_gradient(low = "white", high = "red",
                                      name = "Correlation\n[R or rho]") +
                  xlab("") +
                  ylab("") +
                  labs(title = paste(column_values, "correlation between samples:\nPearson (upper triangle) and Spearman correlation (lower triangle)")) +
                  scale_x_discrete(expand = c(0, 0)) +
                  scale_y_discrete(limits = rev(levels(data.plot$Var1)), expand = c(0, 0)) +
                  theme_bw() +
                  theme(plot.title = element_text(hjust = 0.5, vjust = 1),
                        axis.text.x = element_text(angle = 90, hjust = 1)))
    }
    print(p)
    return(data.plot)
}
