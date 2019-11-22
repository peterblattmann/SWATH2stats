utils::globalVariables(c("sd", "na.omit", "aggregate", "density"))

plot_variation <- function(data, column.values = "Intensity", 
                           Comparison = transition_group_id + Condition ~ BioReplicate, 
                           fun.aggregate = NULL, 
                           label=FALSE, 
                           title = "cv across conditions",
                           boxplot = TRUE, ...){
  if(sum(colnames(data) == "decoy") == 1){
    data <- data[data$decoy == 0,]
  }
  data.c <- dcast(data, Comparison, value.var = column.values, fun.aggregate = fun.aggregate)
  data.c[data.c == 0] <- NA
  
  n_vars <- length(all.vars(Comparison))
  
  data.sd <- apply(data.c[,n_vars:dim(data.c)[2]], 1, function(x) sd(x, na.rm=TRUE))
  data.mean <- apply(data.c[,n_vars:dim(data.c)[2]], 1, function(x) mean(x, na.rm=TRUE))
  data.c$cv <- data.sd/data.mean
  mean.cv <- mean(data.c$cv, na.rm=TRUE)
  median.cv <- median(data.c$cv, na.rm=TRUE)

  data.cv <- data.c[,c(colnames(data.c)[2], "cv")]
  
  p <- ggplot(na.omit(data.cv), aes_string(x=colnames(data.cv)[1], y="cv")) + 
    geom_violin(scale="area") + 
    theme_bw() +
    theme(axis.text.x = element_text(size= 8, angle = 90, hjust = 1, vjust = 0.5)) + 
    labs(title= title)
  
  
  if(isTRUE(label)){
    p <- p + stat_summary(fun.data = function(x)data.frame(y=max(x)*0.75,label=paste("median cv:\n", signif(median(x,na.rm=TRUE), digits=2))), geom="text")

  }
  
  if(isTRUE(boxplot)){
    p <- p + geom_boxplot(width=0.1, outlier.shape = NA)
  }


  print(p)
  if("Condition" %in% colnames(data.c)){
    median <- aggregate(data.c[,"cv"], by=list(data.c$Condition), FUN=function(x)median(x, na.rm=TRUE))
    colnames(median) <- c("Condition", "median_cv")
    mean <- aggregate(data.c[,"cv"], by=list(data.c$Condition), FUN=function(x)mean(x, na.rm=TRUE))
    colnames(mean) <- c("Condition", "mean_cv")
    mode <- aggregate(data.c[,"cv"], by=list(data.c$Condition), FUN=function(x){d<-density(x, na.rm=TRUE); i <- which.max(d$y); return(d$x[i])})
    colnames(mode) <- c("Condition", "mode_cv")

    cv_table <- merge(mode, merge(mean, median, by="Condition"), by="Condition")

    return(list(data.c, cv_table))

  }
  return(data.c)


}
