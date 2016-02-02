utils::globalVariables(c("scope", "cv", "value"))

plot_variation_vs_total <- function(data, column.values = "Intensity", Comparison1 = transition_group_id  ~ BioReplicate + Condition, Comparison2 = transition_group_id + Condition ~ BioReplicate, fun.aggregate = NULL, label=TRUE,...){
    if(sum(colnames(data) == "decoy") == 1){
    data <- data[data$decoy == 0,]
  }
  data1.c <- dcast(data, Comparison1, value.var = column.values, fun.aggregate = fun.aggregate)
  data1.c[data1.c == 0] <- NA
  data1.sd <- apply(data1.c[,2:dim(data1.c)[2]], 1, function(x) sd(x, na.rm=TRUE))
  data1.mean <- apply(data1.c[,2:dim(data1.c)[2]], 1, function(x) mean(x, na.rm=TRUE))
  data1.c$cv <- data1.sd/data1.mean
  mean1.cv <- mean(data1.c$cv, na.rm=TRUE)
  median1.cv <- median(data1.c$cv, na.rm=TRUE)

  data2.c <- dcast(data, Comparison2, value.var = column.values, fun.aggregate = fun.aggregate)
  data2.c[data2.c == 0] <- NA

  data2.sd <- apply(data2.c[,3:dim(data2.c)[2]], 1, function(x) sd(x, na.rm=TRUE))
  data2.mean <- apply(data2.c[,3:dim(data2.c)[2]], 1, function(x) mean(x, na.rm=TRUE))
  data2.c$cv <- data2.sd/data2.mean
  mean2.cv <- mean(data2.c$cv , na.rm=TRUE)
  median2.cv <- median(data2.c$cv , na.rm=TRUE)

  data2.c$rep <- paste(data2.c[,1], data2.c[,2])
  data2.c$scope <- "replicate"
  data1.c$rep <- data1.c[,1]
  data1.c$scope <- "total"

  data.comb <- rbind(data1.c[,c("rep", "cv", "scope")], data2.c[,c("rep", "cv", "scope")])
  data.comb$scope <- factor(data.comb$scope, levels=c("total","replicate"))
  p <- (ggplot(na.omit(data.comb), aes(x=scope, y=cv))
        + geom_violin(scale="area") + xlab("")
        + theme(axis.text.x = element_text(size= 8, angle = 90, hjust = 1, vjust = 0.5))
        + labs(title= "Coefficient of variation - total versus within replicates"))
  if(isTRUE(label)){
    p <- (ggplot(na.omit(data.comb), aes(x=scope, y=cv))
          + geom_violin(scale="area") + xlab("")
          + theme(axis.text.x = element_text(size= 8, angle = 90, hjust = 1, vjust = 0.5))
          + labs(title= "Coefficient of variation - total versus within replicates")
          + stat_summary(fun.data = function(x)data.frame(y=median(x),label=paste("median cv:\n", signif(median(x,na.rm=T), digits=2))), geom="text")
    )
   }
  print(p)
}
