plot_variation <- function(data, column.values = "Intensity", Comparison = transition_group_id + Condition ~ BioReplicate, fun.aggregate = NULL, label=TRUE,...){
  if(sum(colnames(data) == "decoy") == 1){
    data <- data[data$decoy == 0,]
  }
  data.c <- dcast(data, Comparison, value.var = column.values)
  data.c[data.c == 0] <- NA

  data.sd <- apply(data.c[,3:dim(data.c)[2]], 1, function(x) sd(x, na.rm=TRUE))
  data.mean <- apply(data.c[,3:dim(data.c)[2]], 1, function(x) mean(x, na.rm=TRUE))
  data.c$cv <- data.sd/data.mean
  #mean.cv <- mean(data.c$cv, na.rm=TRUE)
  #median.cv <- median(data.c$cv, na.rm=TRUE)

  p <- (ggplot(na.omit(data.c), aes_string(x=colnames(data.c)[2], y="cv"))
        + geom_violin(scale="area")
        + theme(axis.text.x = element_text(size= 8, angle = 90, hjust = 1, vjust = 0.5))
        + labs(title= "CV across conditions"))

  if(isTRUE(label)){
    p <- (ggplot(na.omit(data.c), aes_string(x=colnames(data.c)[2], y="cv"))
          + geom_violin(scale="area")
          + theme(axis.text.x = element_text(size= 8, angle = 90, hjust = 1, vjust = 0.5))
          + stat_summary(fun.data = function(x)data.frame(y=median(x),label=paste("median cv:\n", signif(median(x,na.rm=T), digits=2))), geom="text")
          + labs(title= "CV across conditions"))
  }

  print(p)
}
