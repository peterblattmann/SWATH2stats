#' Plots the total variation versus variation within replicates
#'
#' This function plots the total variation and the variation within replicates
#' for a given value. If decoys are present these are removed before plotting.
#'
#' @param data Data table that is produced by the OpenSWATH/pyProphet workflow.
#' @param column.values  Indicates the columns for which the variation is
#'   assessed. This can be the Intensity or Signal, but also the retention time.
#' @param comparison1 The comparison for assessing the total
#'   variability. Default is to assess the variability per transition_group_id
#'   over the combination of Replicates and different Conditions.
#' @param comparison2  The comparison for assessing the variability within the
#'   replicates. Default is to assess the variability per transition_group_id
#'   and Condition over the different Replicates.
#' @param fun_aggregate If depending on the comparison values have to be
#'   aggregated one needs to provide the function here. (I think this should be
#'   sum, yesno?)
#' @param label  Option to print value of median cv.
#' @param title Title of plot. Default: "cv across conditions"
#' @param boxplot   Logical. If boxplot should be plotted. Default: TRUE
#' @param ...  Arguments passed through, currently unused.
#' @return Plots in Rconsole a violin plot comparing the total variation with
#'   the variation within replicates. In addition it returns the data frame from
#'   which the plotting is done and a table with the calculated mean, median and
#'   mode of the cv for the total or replicate data.
#' @author Peter Blattmann
#' @examples {
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  var_summary <- plot_variation_vs_total(data)
#'  }
#' @importFrom ggplot2 ggplot geom_violin aes_string geom_tile 
#'   scale_fill_gradient xlab ylab labs geom_text geom_boxplot scale_x_discrete 
#'   scale_y_discrete theme element_text stat_summary
#' @importFrom stats sd na.omit aggregate density
#' @importFrom reshape2 dcast
#' @export
plot_variation_vs_total <- function(data, 
                                    column.values = "Intensity", 
                                    comparison1 = transition_group_id ~ BioReplicate + Condition, 
                                    comparison2 = transition_group_id + Condition ~ BioReplicate,
                                    fun_aggregate = NULL, 
                                    label = FALSE, 
                                    title = "coefficient of variation - total versus within replicates",
                                    boxplot = TRUE, ...){
    if (sum(colnames(data) == "decoy") == 1){
        data <- data[data$decoy == 0, ]
    }
    data1.c <- dcast(data, comparison1, value.var = column.values, 
                     fun.aggregate = fun_aggregate)
    data1.c[data1.c == 0] <- NA
    data1.sd <- apply(data1.c[, 2:dim(data1.c)[2]], 1, 
                      function(x) sd(x, na.rm = TRUE))
    data1.mean <- apply(data1.c[, 2:dim(data1.c)[2]], 1, 
                        function(x) mean(x, na.rm = TRUE))
    data1.c$cv <- data1.sd/data1.mean
    mean1.cv <- mean(data1.c$cv, na.rm = TRUE)
    median1.cv <- median(data1.c$cv, na.rm = TRUE)

    data2.c <- dcast(data, comparison2, value.var = column.values, 
                     fun.aggregate = fun_aggregate)
    data2.c[data2.c == 0] <- NA
    n_vars2 <- length(all.vars(comparison2))

    data2.sd <- apply(data2.c[, n_vars2:dim(data2.c)[2]], 1, 
                      function(x) sd(x, na.rm = TRUE))
    data2.mean <- apply(data2.c[, n_vars2:dim(data2.c)[2]], 1, 
                        function(x) mean(x,
        na.rm = TRUE))
    data2.c$cv <- data2.sd/data2.mean
    mean2.cv <- mean(data2.c$cv, na.rm = TRUE)
    median2.cv <- median(data2.c$cv, na.rm = TRUE)

    data2.c$rep <- paste(data2.c[, 1], data2.c[, 2])
    data2.c$scope <- "replicate"
    data1.c$rep <- data1.c[, 1]
    data1.c$scope <- "total"

    data.comb <- rbind(data1.c[, c("rep", "cv", "scope")], data2.c[, c("rep", "cv",
        "scope")])
    data.comb$scope <- factor(data.comb$scope, 
                              levels = c("total", "replicate"))
    p <- (ggplot(na.omit(data.comb), aes_string(x = "scope", y = "cv")) +
          geom_violin(scale = "area") +
          xlab("") +
          theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1,
                                           vjust = 0.5)) +
          labs(title = title))
    if (label) {
        p <- p +
          stat_summary(
            fun.data = function(x) data.frame("y" = max(x) * 0.75,
                                              label = paste("median cv:\n", 
                                                            signif(median(x, na.rm = TRUE), 
                                                                   digits = 2))),
            geom = "text")
    }
    if (boxplot) {
        p <- p + geom_boxplot(width = 0.1, outlier.shape = NA)
    }
    print(p)

    median <- aggregate(data.comb[, "cv"], by = list(data.comb$scope),
                        FUN = function(x) median(x, na.rm = TRUE))
    colnames(median) <- c("scope", "median_cv")
    mean <- aggregate(data.comb[, "cv"], by = list(data.comb$scope),
                      FUN = function(x) mean(x, na.rm = TRUE))
    colnames(mean) <- c("scope", "mean_cv")
    mode <- aggregate(data.comb[, "cv"], by = list(data.comb$scope),
                      FUN = function(x) {
                        d <- density(x, na.rm = TRUE)
                        i <- which.max(d$y)
                        return(d$x[i])
                      })
    colnames(mode) <- c("scope", "mode_cv")

    cv_table <- merge(mode, merge(mean, median, by = "scope"), by = "scope")
    return(list(data.comb, cv_table))
}
