#'  Make a plot of the correlation between samples/groups in the swath2stats data.
#'
#' @param data SWATH2stats data to poke.
#' @param column.values  Which column to query?
#' @param size  How large should the text in the grid be (smaller is better for
#'   rmarkdown html reports).
#' @param comparison  Formula describing the organization of samples to query.
#' @param fun.aggregate  Which function should be used to gather the data?
#' @param label  Put labels of the correlations in the plot?
#' @return plot!
#' @export
plot_correlation_between_samples <- function(data, column.values="intensity", size=6,
                                             comparison=transition_group_id ~ condition + bioreplicate,
                                             fun.aggregate=NULL, label=TRUE, ...) {

  if (sum(colnames(data) == "decoy") == 1) {
    data <- data[data[["decoy"]] == 0, ]
  }

  data.c <- reshape2::dcast(data, comparison, value.var=column.values,
                            fun.aggregate=fun.aggregate)

  dep.vars <- length(all.vars(comparison[[2]]))
  indep.vars <- length(all.vars(comparison[[3]]))

  pearson.cor <- cor(data.c[, (dep.vars + 1):dim(data.c)[2]],
                     use="pairwise.complete.obs",
                     method="pearson")
  pearson.cor[lower.tri(pearson.cor)] <- NA

  spearman.cor <- cor(data.c[, (dep.vars + 1):dim(data.c)[2]],
                      use="pairwise.complete.obs",
                      method="spearman")
  spearman.cor[upper.tri(spearman.cor, diag = TRUE)] <- NA

  pearson.cor <- reshape2::melt(pearson.cor)
  pearson.cor[["method"]] <- "pearson"
  spearman.cor <- reshape2::melt(spearman.cor)
  spearman.cor[["method"]] <- "spearman"

  data.plot <- rbind(pearson.cor, spearman.cor)
  data.plot <- data.plot[!is.na(data.plot$value), ]

  if (isTRUE(label)) {
    p <- ggplot(data.plot, aes(x=Var2, y=Var1, fill=value)) +
      geom_tile() +
      scale_fill_gradient(low="white", high="red", name="Correlation\n[R or rho]") +
      xlab("") + ylab("") +
      labs(title=paste(column.values, "correlation between samples:\nPearson (upper triangle) and Spearman correlation (lower triangle)")) +
      geom_text(size=size, aes(label=round(data.plot$value, digits=2))) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(limits=rev(levels(data.plot$Var1)), expand=c(0,0)) +
      theme(plot.title=element_text(hjust=0.5, vjust=1),
            axis.text.x=element_text(angle=90, hjust=1))
  } else {
    p <- ggplot(data.plot, aes(x=Var2, y=Var1, fill=value)) +
      geom_tile() +
      scale_fill_gradient(low="white", high="red", name="Correlation\n[R or rho]") +
      xlab("") + ylab("") +
      labs(title=paste(column.values, "correlation between samples:\nPearson (upper triangle) and Spearman correlation (lower triangle)")) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(limits=rev(levels(data.plot$Var1)), expand=c(0,0)) +
      theme(plot.title=element_text(hjust=0.5, vjust=1),
            axis.text.x=element_text(angle=90, hjust=1))
  }

  ## I would prefer returning a list containing the data and the plot object,
  ## then leaving it to the user to decide to print the plot, but I am just
  ## weird that way I think.
  print(p)
  return(data.plot)
}
