plot_variation <- function(data, column.values="intensity",
                           comparison=transition_group_id + condition ~ bioreplicate,
                           fun.aggregate=NULL, label=TRUE, ...) {
  if (sum(colnames(data) == "decoy") == 1) {
    data <- data[data[["decoy"]] == 0, ]
  }
  data.c <- reshape2::dcast(data, comparison, value.var=column.values, fun.aggregate=fun.aggregate)
  data.c[data.c == 0] <- NA

  n_vars <- length(all.vars(comparison))

  data.sd <- apply(data.c[, n_vars:dim(data.c)[2]], 1, function(x) sd(x, na.rm=TRUE))
  data.mean <- apply(data.c[, n_vars:dim(data.c)[2]], 1, function(x) mean(x, na.rm=TRUE))
  data.c[["cv"]] <- data.sd / data.mean
  mean.cv <- mean(data.c[["cv"]], na.rm=TRUE)
  median.cv <- median(data.c[["cv"]], na.rm=TRUE)

  data.cv <- data.c[, c(colnames(data.c)[2], "cv")]

  if (isTRUE(label)) {
    p <- (ggplot(na.omit(data.cv), aes_string(x=colnames(data.cv)[1], y="cv")) +
          geom_violin(scale="area") +
          theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5)) +
          stat_summary(fun.data=function(x) {
            data.frame(y=median(x),
                       label=paste("median cv:\n",
                                   signif(median(x,na.rm=TRUE),
                                          digits=2))) }, geom="text") +
          labs(title=paste(column.values, "cv across conditions")))
  } else {
    p <- (ggplot(na.omit(data.cv), aes_string(x=colnames(data.cv)[1], y="cv")) +
          geom_violin(scale="area") +
          theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5)) +
          labs(title="cv across conditions"))
  }

  print(p)
  if ("Condition" %in% colnames(data.c)) {
    median <- aggregate(data.c[, "cv"], by=list(data.c[["condition"]]),
                        FUN=function(x) median(x, na.rm=TRUE))
    colnames(median) <- c("condition", "median_cv")
    mean <- aggregate(data.c[, "cv"], by=list(data.c[["condition"]]),
                      FUN=function(x) mean(x, na.rm=TRUE))
    colnames(mean) <- c("condition", "mean_cv")
    mode <- aggregate(data.c[, "cv"], by=list(data.c[["condition"]]),
                      FUN=function(x) {
                        d <- density(x, na.rm=TRUE); i <- which.max(d$y); return(d$x[i]) })
    colnames(mode) <- c("condition", "mode_cv")

    cv_table <- merge(mode, merge(mean, median, by="condition"), by="condition")

    return(list(data.c, cv_table))
  }
  return(data.c)
}
