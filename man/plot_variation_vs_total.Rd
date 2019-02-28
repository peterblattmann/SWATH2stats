\name{plot_variation_vs_total}
\alias{plot_variation_vs_total}
\title{
Plots the total variation versus variation within replicates
}
\description{
This function plots the total variation and the variation within replicates for a given value. If decoys are present these are removed before plotting.
}
\usage{
plot_variation_vs_total(data, column.values = "Intensity", 
Comparison1 = transition_group_id  ~ BioReplicate + Condition, 
Comparison2 = transition_group_id + Condition ~ BioReplicate, 
fun.aggregate = NULL, label=FALSE, 
title = "coefficient of variation - total versus within replicates",
boxplot = TRUE,...)
}
\arguments{
  \item{data}{
Data table that is produced by the OpenSWATH/pyProphet workflow
}
  \item{column.values}{
Indicates the columns for which the variation is assessed. This can be the Intensity or Signal, but also the retention time.
}
   \item{Comparison1}{
The comparison for assessing the total variability. Default is to assess the variability per transition_group_id over the combination of Replicates and different Conditions.}
  \item{Comparison2}{
The comparison for assessing the variability within the replicates. Default is to assess the variability per transition_group_id and Condition over the different Replicates.}
  \item{fun.aggregate}{
If depending on the comparison values have to be aggregated one needs to provide the function here.}
  \item{label}{
Option to print value of median cv.}
  \item{title}{
Option to define title of plot.}
  \item{boxplot}{
Option to add boxplot.}
  \item{\dots}{
further arguments passed to method.
}
}
\value{
Plots in Rconsole a violin plot comparing the total variation with the variation within replicates. In addition it returns the data frame from which the plotting is done and a table with the calculated mean, median and mode of the cv for the total or replicate data.
}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
plot_variation_vs_total(data)
}
