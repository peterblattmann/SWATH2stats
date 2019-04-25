utils::globalVariables(c("Var1", "Var2", "value"))

#'  Plots the correlation between injections.
#'
#' This function plots the Pearson's and Spearman correlation between
#' samples. If decoys are present these are removed before plotting.
#'
#' @param data Data frame that is produced by the OpenSWATH/pyProphet workflow.
#' @param column.values   Indicates the columns for which the correlation is
#'   assessed. This can be the Intensity or Signal, but also the retention
#'   time.
#' @param size  How large should the text in the grid be (smaller is better for
#'   rmarkdown html reports).
#' @param comparison  The comparison for assessing the variability. Default is
#'   to assess the variability per transition_group_id over the different
#'   Condition and Replicates. Comparison is performed using the dcast()
#'   function of the reshape2 package.
#' @param fun.aggregate  If for the comparison values have to be aggregated one
#'   needs to provide the function here.
#' @param label  Option to print correlation value in the plot.
#' @param ... Further arguments passed to methods.
#' @return Plots in Rconsole a correlation heatmap and returns the data frame
#'   used to do the plotting.
#' @author Peter Blattmann
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  information <- plot_correlation_between_samples(data)
#' @export
plot_correlation_between_samples <- function(data, column.values = "Intensity",
                                             Comparison = transition_group_id ~ Condition + BioReplicate,
                                             fun.aggregate = NULL, label = TRUE, ...) {
    if (sum(colnames(data) == "decoy") == 1) {
        data <- data[data$decoy == 0, ]
    }
    data.c <- dcast(data, Comparison, value.var = column.values, fun.aggregate = fun.aggregate)
    dep.vars <- length(all.vars(Comparison[[2]]))
    indep.vars <- length(all.vars(Comparison[[3]]))

    pearson.cor <- cor(data.c[, (dep.vars + 1):dim(data.c)[2]], use = "pairwise.complete.obs",
        method = "pearson")
    pearson.cor[lower.tri(pearson.cor)] <- NA

    spearman.cor <- cor(data.c[, (dep.vars + 1):dim(data.c)[2]], use = "pairwise.complete.obs",
        method = "spearman")
    spearman.cor[upper.tri(spearman.cor, diag = TRUE)] <- NA

    pearson.cor <- melt(pearson.cor)
    pearson.cor$method <- "pearson"
    spearman.cor <- melt(spearman.cor)
    spearman.cor$method <- "spearman"

    data.plot <- rbind(pearson.cor, spearman.cor)
    data.plot <- data.plot[!is.na(data.plot$value), ]

    if (isTRUE(label)) {
        p <- (ggplot(data.plot, aes_string(x = "Var2", y = "Var1", fill = "value")) +
              ggplot2::geom_tile() +
              ggplot2::scale_fill_gradient(low = "white", high = "red",
                                           name = "Correlation\n[R or rho]") +
              ggplot2::xlab("") +
              ggplot2::ylab("") +
              ggplot2::labs(title = paste(column.values, "correlation between samples:\nPearson (upper triangle) and Spearman correlation (lower triangle)")) +
              ggplot2::geom_text(aes_string(label = "round(data.plot$value, digits = 2)")) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_discrete(limits = rev(levels(data.plot$Var1)),
                                        expand = c(0, 0)) +
              ggplot2::theme(plot.title = element_text(hjust = 0.5, vjust = 1),
                             axis.text = ggplot2::element_text(angle = 90, hjust = 1)))
    } else {
        p <- (ggplot(data.plot, aes_string(x = "Var2", y = "Var1", fill = "value")) +
              ggplot2::geom_tile() +
              ggplot2::scale_fill_gradient(low = "white", high = "red",
                                           name = "Correlation\n[R or rho]") +
              ggplot2::xlab("") +
              ggplot2::ylab("") +
              ggplot2::labs(title = paste(column.values, "correlation between samples:\nPearson (upper triangle) and Spearman correlation (lower triangle)")) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_discrete(limits = rev(levels(data.plot$Var1)), expand = c(0, 0)) +
              ggplot2::theme(plot.title = element_text(hjust = 0.5, vjust = 1),
                             axis.text = ggplot2::element_text(angle = 90, hjust = 1)))
    }
    print(p)
    return(data.plot)
}
