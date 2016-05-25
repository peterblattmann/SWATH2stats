\name{plot_correlation_between_samples}
\alias{plot_correlation_between_samples}
\title{
Plots the correlation between injections.
}
\description{
This function plots the Pearson's and Spearman correlation between samples. If decoys are present these are removed before plotting.
}
\usage{
plot_correlation_between_samples(data, column.values = "Intensity", 
Comparison = transition_group_id ~ Condition + BioReplicate, 
fun.aggregate =NULL, label=TRUE, ...)
}
\arguments{
  \item{data}{
Data frame that is produced by the OpenSWATH/pyProphet workflow
}
  \item{column.values}{
Indicates the columns for which the correlation is assessed. This can be the Intensity or Signal, but also the retention time.
}
  \item{Comparison}{
The comparison for assessing the variability. Default is to assess the variability per transition_group_id over the different Condition and Replicates. Comparison is performed using the dcast() function of the reshape2 package.}
    \item{fun.aggregate}{
If for the comparison values have to be aggregated one needs to provide the function here.}
    \item{label}{
Option to print correlation value in the plot.}
  \item{\dots}{
further arguments passed to method.}
}
\value{
Plots in Rconsole a correlation heatmap and returns the data frame used to do the plotting.
}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
plot_correlation_between_samples(data)
}
