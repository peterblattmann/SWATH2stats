utils::globalVariables(c("Var1", "Var2", "value"))

plot_correlation_between_samples <- function(data, column.values = "Intensity", Comparison = transition_group_id ~ Condition + BioReplicate, ...){

  if(sum(colnames(data) == "decoy") == 1){
    data <- data[data$decoy == 0,]
  }

  data.c <- dcast(data, Comparison, value.var = column.values, fun.aggregate=sum)

  pearson.cor <- cor(data.c[,2:dim(data.c)[2]], use="pairwise.complete.obs", method="pearson")
  pearson.cor[lower.tri(pearson.cor)] <- NA

  spearman.cor <- cor(data.c[,2:dim(data.c)[2]], use="pairwise.complete.obs", method="spearman")
  spearman.cor[upper.tri(spearman.cor, diag = TRUE)] <- NA


  data.plot <- rbind(melt(pearson.cor), melt(spearman.cor))
  data.plot <- data.plot[!is.na(data.plot$value),]

  p <- (ggplot(data.plot, aes(x=Var2, y=Var1, fill=value)) + geom_tile()
        + scale_fill_gradient(low = "white", high="red", name="Correlation\n[R or rho]")
        + xlab("") + ylab("")
        + labs(title="Correlation between samples: Pearson  (upper triangle) and Spearman correlation (lower triangle)")
        + geom_text(aes(fill = data.plot$value, label = round(data.plot$value, digits= 2)))
        + theme(plot.title = element_text(hjust = 0, vjust = 1)))

  print(p)

}
